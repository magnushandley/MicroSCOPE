import argparse
import sys

import uproot
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import shap

from typing import Dict, List, Optional

from sklearn.metrics import roc_auc_score, classification_report

import xgboost as xgb
from pathlib import Path

SCORE_COLUMN = "bdt_score"


def source_files_in_order(df: pd.DataFrame) -> List[str]:
    """Return source files in their first-seen row order."""
    if "source_file" not in df.columns:
        raise ValueError("DataFrame is missing required metadata column 'source_file'")
    return df["source_file"].drop_duplicates().astype(str).tolist()


def default_train_fractions(root_files: List[str], test_size: float) -> List[float]:
    """Return the legacy uniform train fraction derived from --test-size."""
    return [1.0 - test_size] * len(root_files)


def validate_train_fractions(
    fractions: List[float],
    root_files: List[str],
    arg_name: str,
    file_arg_name: str,
) -> None:
    if len(fractions) != len(root_files):
        raise ValueError(
            f"{arg_name} must have the same length as {file_arg_name}: "
            f"{len(fractions)} vs {len(root_files)}"
        )

    invalid = [
        fraction
        for fraction in fractions
        if not 0.0 <= fraction < 1.0
    ]
    if invalid:
        raise ValueError(f"{arg_name} values must be >= 0.0 and < 1.0")


def normalize_data_files(
    data_file: Optional[str],
    data_files: Optional[List[str]],
) -> List[str]:
    """Merge legacy --data-file with plural --data-files, preserving order."""
    normalized = []
    seen = set()

    for fname in data_files or []:
        if fname not in seen:
            normalized.append(fname)
            seen.add(fname)

    if data_file and data_file not in seen:
        normalized.append(data_file)

    return normalized


def train_fraction_map(
    root_files: List[str],
    train_fractions: List[float],
) -> Dict[str, float]:
    return {
        root_file: train_fraction
        for root_file, train_fraction in zip(root_files, train_fractions)
    }


def load_root_files_to_df(
    root_files: List[str],
    tree_name: str,
    branches: Optional[List[str]] = None,
    library: str = "pd",
    per_file_weights: Optional[List[float]] = None,
    weight_column: str = "sample_weight",
) -> pd.DataFrame:
    """
    Load multiple ROOT files into a single pandas DataFrame.

    Parameters
    ----------
    root_files : list[str]
        List of paths to ROOT files.
    tree_name : str
        Name of the TTree (e.g. "myTree" or "events").
    branches : list[str] or None
        Branches to read. If None, read all branches.
    library : str
        Backend for uproot ("pd" for pandas, "np" for numpy).
    per_file_weights : list[float] or None
        Optional list of weights, one per ROOT file, applied to every row read from that file.
    weight_column : str
        Name of the column to store the per-row weight.

    Returns
    -------
    pandas.DataFrame
    """
    dfs = []

    if per_file_weights is not None and len(per_file_weights) != len(root_files):
        raise ValueError(
            f"per_file_weights must have same length as root_files: "
            f"{len(per_file_weights)} vs {len(root_files)}"
        )

    for i, fname in enumerate(root_files):
        print(f"Loading {fname}")
        with uproot.open(fname) as f:
            tree = f[tree_name]

            df = tree.arrays(
                expressions=branches,
                library=library
            )

            # Stable identifier for matching rows back to the original file
            df["row_in_file"] = np.arange(len(df), dtype=np.int64)

            if per_file_weights is not None:
                df[weight_column] = float(per_file_weights[i])
            else:
                df[weight_column] = 1.0

            df["source_file"] = fname

            dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def prepare_train_test(
    df_bkg: pd.DataFrame,
    df_sig: pd.DataFrame,
    feature_columns: List[str],
    weight_column: str = "sample_weight",
    bkg_train_fractions: Optional[List[float]] = None,
    sig_train_fractions: Optional[List[float]] = None,
    test_size: float = 0.5,
    random_state: int = 1337,
    balance_sig_to_bkg: bool = True,
):
    """Prepare X/y/weights and split into train/test per input file.

    Weighting strategy:
      - Background component imbalance is handled by the per-row weights already assigned
        from `bkg_sample_weights` (one per background file).
      - Optionally balance the *total weighted* signal yield to the *total weighted* background yield
        by scaling signal weights so sum_w_sig == sum_w_bkg (overall factor).
    """

    # Labels
    df_bkg = df_bkg.copy()
    df_sig = df_sig.copy()
    df_bkg["label"] = 0
    df_sig["label"] = 1

    # Ensure weight column exists
    if weight_column not in df_bkg.columns:
        df_bkg[weight_column] = 1.0
    if weight_column not in df_sig.columns:
        df_sig[weight_column] = 1.0

    bkg_files = source_files_in_order(df_bkg)
    sig_files = source_files_in_order(df_sig)
    if bkg_train_fractions is None:
        bkg_train_fractions = default_train_fractions(bkg_files, test_size)
    if sig_train_fractions is None:
        sig_train_fractions = default_train_fractions(sig_files, test_size)

    validate_train_fractions(
        bkg_train_fractions,
        bkg_files,
        "--bkg-train-fractions",
        "--bkg-files",
    )
    validate_train_fractions(
        sig_train_fractions,
        sig_files,
        "--sig-train-fractions",
        "--sig-files",
    )

    # Optional: scale signal weights to match total weighted background yield
    if balance_sig_to_bkg:
        sum_w_bkg = float(df_bkg[weight_column].sum())
        sum_w_sig = float(df_sig[weight_column].sum())
        if sum_w_sig > 0:
            scale = sum_w_bkg / sum_w_sig
            df_sig[weight_column] = df_sig[weight_column] * scale
            print(f"[weights] Scaling signal weights by {scale:.6g} so sum_w_sig == sum_w_bkg")
        else:
            print("[weights] Warning: sum_w_sig is 0; not scaling signal weights")

    # Combine
    df = pd.concat([df_bkg, df_sig], ignore_index=True)

    train_fractions = {}
    train_fractions.update(train_fraction_map(bkg_files, bkg_train_fractions))
    train_fractions.update(train_fraction_map(sig_files, sig_train_fractions))

    train_indices = []
    test_indices = []
    ordered_files = bkg_files + sig_files
    for file_index, source_file in enumerate(ordered_files):
        file_indices = df.index[df["source_file"] == source_file].to_numpy()
        fraction = train_fractions[source_file]
        rng = np.random.default_rng(random_state + file_index)
        shuffled_indices = rng.permutation(file_indices)
        n_train = int(np.floor(fraction * len(shuffled_indices)))

        train_indices.extend(shuffled_indices[:n_train].tolist())
        test_indices.extend(shuffled_indices[n_train:].tolist())

    # Build arrays
    X = df[feature_columns]
    y = df["label"].astype(int)
    w = df[weight_column].astype(float)

    X_train = X.loc[train_indices]
    X_test = X.loc[test_indices]
    y_train = y.loc[train_indices]
    y_test = y.loc[test_indices]
    w_train = w.loc[train_indices]
    w_test = w.loc[test_indices]

    for label, class_name in [(0, "background"), (1, "signal")]:
        if not (y_train == label).any():
            raise ValueError(f"Training split contains no {class_name} rows")
        if not (y_test == label).any():
            raise ValueError(f"Held-out split contains no {class_name} rows")

    # Keep metadata to allow filtering back to original files
    meta_cols = ["source_file", "row_in_file", "label", weight_column]
    meta = df[meta_cols]
    meta_train = meta.loc[train_indices].copy()
    meta_test = meta.loc[test_indices].copy()

    return X_train, X_test, y_train, y_test, w_train, w_test, meta_train, meta_test


def train_xgboost_bdt(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    w_train: pd.Series,
    X_test: pd.DataFrame,
    y_test: pd.Series,
    w_test: pd.Series,
    *,
    random_state: int = 1337,
    num_boost_round: int = 5000,
    early_stopping_rounds: int = 50,
):
    """Train an XGBoost BDT using native xgb.train() for broad version compatibility.

    This avoids sklearn-wrapper API differences (e.g. missing early_stopping_rounds/callbacks).
    Per-event weights are applied via DMatrix weights.
    """

    # Convert to DMatrix (native XGBoost format)
    dtrain = xgb.DMatrix(X_train, label=y_train, weight=w_train)
    dtest = xgb.DMatrix(X_test, label=y_test, weight=w_test)

    params = {
        "objective": "binary:logistic",
        "eval_metric": "auc",
        "max_depth": 4,
        "eta": 0.05,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "min_child_weight": 1.0,
        "lambda": 1.0,
        "gamma": 0.0,
        "seed": random_state,
        "verbosity": 1,
        # Change to "gpu_hist" if your installation supports GPU
        "tree_method": "hist",
    }

    evals = [(dtrain, "train"), (dtest, "eval")]

    booster = xgb.train(
        params=params,
        dtrain=dtrain,
        num_boost_round=num_boost_round,
        evals=evals,
        early_stopping_rounds=early_stopping_rounds,
        verbose_eval=False,
    )

    # Best iteration info
    best_iter = getattr(booster, "best_iteration", None)
    best_score = getattr(booster, "best_score", None)
    if best_iter is not None:
        print(f"[xgboost] Best iteration: {best_iter}")
    if best_score is not None:
        print(f"[xgboost] Best score: {best_score}")

    # Predictions on test
    y_proba = booster.predict(dtest)
    auc = roc_auc_score(y_test, y_proba, sample_weight=w_test)
    print(f"[xgboost] Weighted ROC AUC (test): {auc:.6f}")

    y_pred = (y_proba >= 0.5).astype(int)
    print("[xgboost] Classification report (unweighted, threshold=0.5):")
    print(classification_report(y_test, y_pred, digits=4))

    return booster


def print_shap_importance_table(
    shap_values,
    feature_names: List[str],
    threshold: float = 0.999,
) -> None:
    """Print ranked mean absolute SHAP importance and a cumulative cutoff."""
    values = np.asarray(shap_values.values)
    if values.ndim != 2:
        raise ValueError("SHAP values must have shape (events, features)")
    if values.shape[1] != len(feature_names):
        raise ValueError(
            "SHAP feature count does not match the supplied feature names: "
            f"{values.shape[1]} vs {len(feature_names)}"
        )
    if not 0.0 < threshold <= 1.0:
        raise ValueError("SHAP importance threshold must be in (0, 1]")

    mean_abs_shap = np.mean(np.abs(values), axis=0)
    order = np.argsort(-mean_abs_shap, kind="stable")
    sorted_importance = mean_abs_shap[order]
    total_importance = float(sorted_importance.sum())

    keep = np.zeros(len(feature_names), dtype=bool)
    if total_importance > 0.0:
        fractional_importance = sorted_importance / total_importance
        cumulative_importance = np.cumsum(fractional_importance)
        n_keep = int(np.searchsorted(cumulative_importance, threshold)) + 1
        keep[:n_keep] = True
    else:
        fractional_importance = np.zeros(len(feature_names), dtype=float)
        cumulative_importance = np.zeros(len(feature_names), dtype=float)

    keep_column = f"Keep at {threshold:.1%}"
    table = pd.DataFrame(
        {
            "Rank": np.arange(1, len(feature_names) + 1),
            "Variable": np.asarray(feature_names)[order],
            "Mean |SHAP|": sorted_importance,
            "Fractional importance (%)": 100.0 * fractional_importance,
            "Cumulative importance (%)": 100.0 * cumulative_importance,
            keep_column: keep,
        }
    )

    print(f"[xgboost] SHAP feature importance ({threshold:.1%} cutoff):")
    print(
        table.to_string(
            index=False,
            formatters={
                "Mean |SHAP|": "{:.6g}".format,
                "Fractional importance (%)": "{:.6f}".format,
                "Cumulative importance (%)": "{:.6f}".format,
            },
        )
    )
    if total_importance == 0.0:
        print(
            "[xgboost] WARNING: all mean absolute SHAP values are zero; "
            "no meaningful cumulative cutoff can be calculated."
        )


def save_shap_beeswarm(
    booster: xgb.Booster,
    X_test: pd.DataFrame,
    output_path: Path,
    random_state: int,
    max_events: int = 10_000,
) -> None:
    """Save a SHAP beeswarm diagnostic for a deterministic test-set sample."""
    if X_test.empty:
        raise ValueError("Cannot make a SHAP beeswarm from an empty test set")

    if len(X_test) > max_events:
        X_shap = X_test.sample(n=max_events, random_state=random_state)
    else:
        X_shap = X_test.copy()

    explainer = shap.TreeExplainer(booster)
    shap_values = explainer(X_shap)
    print_shap_importance_table(shap_values, list(X_shap.columns))
    shap.plots.beeswarm(
        shap_values,
        max_display=X_shap.shape[1],
        show=False,
        plot_size="auto",
    )

    figure = plt.gcf()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"[xgboost] SHAP beeswarm events: {len(X_shap)}")
    print(f"[xgboost] Saved SHAP beeswarm to {output_path}")


# --- Helper functions for scoring and writing ROOT ---
def add_bdt_score(
    booster: xgb.Booster,
    df: pd.DataFrame,
    feature_columns: List[str],
    score_column: str = SCORE_COLUMN,
) -> pd.DataFrame:
    """Return a copy of df with an added BDT score column."""
    df_out = df.copy()
    dmat = xgb.DMatrix(df_out[feature_columns])
    df_out[score_column] = booster.predict(dmat)
    return df_out


def dedupe_preserve_order(branches: List[str]) -> List[str]:
    """Return branch names without duplicates, preserving the first occurrence."""
    seen = set()
    unique = []
    for branch in branches:
        if branch in seen:
            continue
        seen.add(branch)
        unique.append(branch)
    return unique


def validate_no_reserved_score_branch(branches: List[str], arg_name: str) -> None:
    if SCORE_COLUMN in branches:
        raise ValueError(
            f"{arg_name} cannot contain reserved output branch '{SCORE_COLUMN}'"
        )


def resolve_output_branches_for_tree(
    tree,
    required_branches: List[str],
    optional_branches: Optional[List[str]],
    in_path: str,
) -> List[str]:
    """Return output branches present in this tree, failing only for required branches."""
    available = set(tree.keys())
    required_unique = dedupe_preserve_order(required_branches)
    optional_unique = [
        branch
        for branch in dedupe_preserve_order(optional_branches or [])
        if branch not in required_unique
    ]

    missing_required = [
        branch for branch in required_unique
        if branch not in available
    ]
    if missing_required:
        raise ValueError(
            f"Missing required feature branch(es) in {in_path}: "
            f"{', '.join(missing_required)}"
        )

    present_optional = [
        branch for branch in optional_unique
        if branch in available
    ]
    missing_optional = [
        branch for branch in optional_unique
        if branch not in available
    ]
    if missing_optional:
        print(
            "\n"
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
            "[xgboost] WARNING: requested optional output branches are missing\n"
            f"[xgboost] File: {in_path}\n"
            f"[xgboost] Skipping branch(es): {', '.join(missing_optional)}\n"
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
            file=sys.stderr,
        )

    return required_unique + present_optional


def filter_branch_arrays(branch_arrays: dict, allowed_rows: Optional[set]) -> dict:
    """Apply the held-out-row selection to numpy/awkward arrays by entry index."""
    if allowed_rows is None:
        return branch_arrays

    n_entries = len(next(iter(branch_arrays.values()))) if branch_arrays else 0
    row_mask = np.fromiter(
        (row in allowed_rows for row in range(n_entries)),
        dtype=bool,
        count=n_entries,
    )
    return {
        branch_name: branch_array[row_mask]
        for branch_name, branch_array in branch_arrays.items()
    }


def load_output_branch_arrays(
    in_path: str,
    tree_name: str,
    required_branches: List[str],
    optional_branches: Optional[List[str]] = None,
    allowed_rows: Optional[set] = None,
) -> dict:
    """Read branch arrays for output preservation using awkward-compatible uproot arrays."""
    with uproot.open(in_path) as f:
        tree = f[tree_name]
        branches = resolve_output_branches_for_tree(
            tree,
            required_branches,
            optional_branches,
            in_path,
        )
        branch_arrays = tree.arrays(
            expressions=branches,
            library="ak",
            how=dict,
        )

    return filter_branch_arrays(branch_arrays, allowed_rows)


def write_trees_to_root(
    out_path: str,
    trees: dict,
):
    """Write multiple TTrees into one ROOT file.

    Parameters
    ----------
    out_path : str
        Output ROOT filename.
    trees : dict
        Mapping of tree_name -> dict of branch name to numpy/awkward array.
    """

    with uproot.recreate(out_path) as f:
        for tree_name, branch_arrays in trees.items():
            tree_path = [part for part in tree_name.split("/") if part]
            tree_dir = f
            for directory in tree_path[:-1]:
                tree_dir = tree_dir.mkdir(directory)

            tree_leaf_name = tree_path[-1]
            branch_types = {
                branch_name: infer_mktree_branch_type(branch_array)
                for branch_name, branch_array in branch_arrays.items()
            }
            writable_tree = tree_dir.mktree(tree_leaf_name, branch_types)
            writable_tree.extend(branch_arrays)


def infer_mktree_branch_type(branch_array):
    """Infer an uproot mktree branch type from an in-memory branch array."""
    if isinstance(branch_array, np.ndarray):
        if branch_array.ndim <= 1:
            return branch_array.dtype
        fixed_shape = " * ".join(str(dim) for dim in branch_array.shape[1:])
        return f"{fixed_shape} * {branch_array.dtype.name}"

    array_type = getattr(branch_array, "type", None)
    content_type = getattr(array_type, "content", None)
    if content_type is not None:
        return str(content_type)
    if array_type is not None:
        return str(array_type)

    return np.asarray(branch_array).dtype


def score_single_root_file(
    booster: xgb.Booster,
    in_path: str,
    out_path: str,
    tree_name: str,
    feature_columns: List[str],
    keep_branches: Optional[List[str]] = None,
    allowed_rows: Optional[set] = None,
):
    """Load one ROOT file, add bdt_score branch, and write to a new ROOT file.

    The output file will contain a TTree with the same name as `tree_name`.
    Requested output branches are read separately from the input file so scalar
    and vector branches can be preserved without routing them through pandas.
    """
    df = load_root_files_to_df(
        [in_path],
        tree_name,
        feature_columns,
        per_file_weights=None,
        weight_column="sample_weight",
    )

    if allowed_rows is not None:
        # Keep only events that were part of the held-out test set
        before = len(df)
        df = df[df["row_in_file"].isin(allowed_rows)].copy()
        after = len(df)
        print(f"[xgboost] Filtering {in_path} to test entries: {before} -> {after}")

    df_scored = add_bdt_score(
        booster,
        df,
        feature_columns=feature_columns,
        score_column=SCORE_COLUMN,
    )

    branch_arrays = load_output_branch_arrays(
        in_path,
        tree_name,
        required_branches=feature_columns,
        optional_branches=keep_branches,
        allowed_rows=allowed_rows,
    )
    branch_arrays[SCORE_COLUMN] = df_scored[SCORE_COLUMN].to_numpy()

    n_scores = len(branch_arrays[SCORE_COLUMN])
    mismatched = [
        branch_name
        for branch_name, branch_array in branch_arrays.items()
        if len(branch_array) != n_scores
    ]
    if mismatched:
        raise ValueError(
            f"Output branch length mismatch for {in_path}: "
            f"{', '.join(mismatched)} do not match {SCORE_COLUMN}"
        )

    # Write a single tree with the original tree_name
    write_trees_to_root(out_path, {tree_name: branch_arrays})

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Train an XGBoost BDT and write scored ROOT files."
    )
    parser.add_argument("--config", help="Optional MicroSCOPE config path passed through by PythonModule.")
    parser.add_argument("--bkg-files", nargs="+", required=True, help="Background ROOT input files.")
    parser.add_argument("--sig-files", nargs="+", required=True, help="Signal ROOT input files.")
    parser.add_argument("--data-file", help="Legacy single data ROOT input file to score.")
    parser.add_argument("--data-files", nargs="+", help="Data ROOT input files to score.")
    parser.add_argument(
        "--det-var-files",
        "--detvar-files",
        nargs="*",
        default=[],
        help="Detector-variation ROOT input files to score like data.",
    )
    parser.add_argument(
        "--bkg-sample-weights",
        nargs="+",
        type=float,
        required=True,
        help="Per-file weights for background files.",
    )
    parser.add_argument(
        "--sig-sample-weights",
        nargs="+",
        type=float,
        required=True,
        help="Per-file weights for signal files.",
    )
    parser.add_argument("--tree-name", required=True, help="Input and output TTree name.")
    parser.add_argument("--branches", nargs="+", required=True, help="Feature branches to train on.")
    parser.add_argument(
        "--keep-branches",
        nargs="*",
        default=[],
        help=(
            "Additional branches to preserve in scored ROOT outputs. "
            "These are written alongside --branches and bdt_score, but are not used for training."
        ),
    )
    parser.add_argument("--output-dir", default=".", help="Directory for model and scored ROOT outputs.")
    parser.add_argument("--model-output", default="xgb_bdt.json", help="Model filename or path.")
    parser.add_argument("--test-size", type=float, default=0.5, help="Held-out test fraction.")
    parser.add_argument(
        "--bkg-train-fractions",
        nargs="+",
        type=float,
        help="Per-background-file training fractions. Defaults to 1 - --test-size.",
    )
    parser.add_argument(
        "--sig-train-fractions",
        nargs="+",
        type=float,
        help="Per-signal-file training fractions. Defaults to 1 - --test-size.",
    )
    parser.add_argument("--random-state", type=int, default=1337, help="Random seed.")
    parser.add_argument("--num-boost-round", type=int, default=5000, help="Maximum XGBoost boosting rounds.")
    parser.add_argument(
        "--early-stopping-rounds",
        type=int,
        default=50,
        help="XGBoost early stopping rounds.",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.bkg_files:
        raise ValueError("--bkg-files must contain at least one file")
    if not args.sig_files:
        raise ValueError("--sig-files must contain at least one file")
    if not normalize_data_files(args.data_file, args.data_files):
        raise ValueError("at least one of --data-file or --data-files must be provided")
    if not args.tree_name:
        raise ValueError("--tree-name cannot be empty")
    if not args.branches:
        raise ValueError("--branches must contain at least one branch")
    validate_no_reserved_score_branch(args.branches, "--branches")
    validate_no_reserved_score_branch(args.keep_branches, "--keep-branches")
    if len(args.bkg_sample_weights) != len(args.bkg_files):
        raise ValueError(
            "--bkg-sample-weights must have the same length as --bkg-files: "
            f"{len(args.bkg_sample_weights)} vs {len(args.bkg_files)}"
        )
    if len(args.sig_sample_weights) != len(args.sig_files):
        raise ValueError(
            "--sig-sample-weights must have the same length as --sig-files: "
            f"{len(args.sig_sample_weights)} vs {len(args.sig_files)}"
        )
    if not 0.0 < args.test_size < 1.0:
        raise ValueError("--test-size must be between 0 and 1")
    if args.bkg_train_fractions is not None:
        validate_train_fractions(
            args.bkg_train_fractions,
            args.bkg_files,
            "--bkg-train-fractions",
            "--bkg-files",
        )
    if args.sig_train_fractions is not None:
        validate_train_fractions(
            args.sig_train_fractions,
            args.sig_files,
            "--sig-train-fractions",
            "--sig-files",
        )
    if args.num_boost_round <= 0:
        raise ValueError("--num-boost-round must be positive")
    if args.early_stopping_rounds <= 0:
        raise ValueError("--early-stopping-rounds must be positive")


def output_path(output_dir: Path, requested_path: str) -> Path:
    path = Path(requested_path)
    if path.is_absolute():
        return path
    return output_dir / path


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        validate_args(args)
    except ValueError as exc:
        parser.error(str(exc))

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_files = normalize_data_files(args.data_file, args.data_files)
    bkg_train_fractions = (
        args.bkg_train_fractions
        if args.bkg_train_fractions is not None
        else default_train_fractions(args.bkg_files, args.test_size)
    )
    sig_train_fractions = (
        args.sig_train_fractions
        if args.sig_train_fractions is not None
        else default_train_fractions(args.sig_files, args.test_size)
    )

    print("[xgboost] Required output feature branches:")
    for branch in dedupe_preserve_order(args.branches):
        print(f"  {branch}")
    if args.keep_branches:
        print("[xgboost] Optional output keep branches:")
        for branch in dedupe_preserve_order(args.keep_branches):
            print(f"  {branch}")
    print("[xgboost] Always-added output branch:")
    print(f"  {SCORE_COLUMN}")
    print("[xgboost] Background train fractions:")
    for fname, fraction in zip(args.bkg_files, bkg_train_fractions):
        print(f"  {fname}: {fraction:.6g}")
    print("[xgboost] Signal train fractions:")
    for fname, fraction in zip(args.sig_files, sig_train_fractions):
        print(f"  {fname}: {fraction:.6g}")

    df_bkg = load_root_files_to_df(
        args.bkg_files,
        args.tree_name,
        args.branches,
        per_file_weights=args.bkg_sample_weights,
        weight_column="sample_weight",
    )

    df_sig = load_root_files_to_df(
        args.sig_files,
        args.tree_name,
        args.branches,
        per_file_weights=args.sig_sample_weights,
        weight_column="sample_weight",
    )

    print(df_bkg.head())
    print(f"Loaded background rows: {len(df_bkg)}")
    print(df_sig.head())
    print(f"Loaded signal rows: {len(df_sig)}")

    X_train, X_test, y_train, y_test, w_train, w_test, meta_train, meta_test = prepare_train_test(
        df_bkg,
        df_sig,
        feature_columns=args.branches,
        weight_column="sample_weight",
        bkg_train_fractions=bkg_train_fractions,
        sig_train_fractions=sig_train_fractions,
        test_size=args.test_size,
        random_state=args.random_state,
        balance_sig_to_bkg=True,
    )

    print(f"Train: n={len(X_train)}  (sum_w={w_train.sum():.6g})")
    print(f"Test:  n={len(X_test)}  (sum_w={w_test.sum():.6g})")

    model = train_xgboost_bdt(
        X_train,
        y_train,
        w_train,
        X_test,
        y_test,
        w_test,
        random_state=args.random_state,
        num_boost_round=args.num_boost_round,
        early_stopping_rounds=args.early_stopping_rounds,
    )

    model_path = output_path(out_dir, args.model_output)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    model.save_model(str(model_path))
    print(f"[xgboost] Saved model to {model_path}")

    #save_shap_beeswarm(
    #    model,
    #    X_test,
    #    out_dir / "shap_beeswarm.png",
    #    random_state=args.random_state,
    #)

    test_rows_by_file = {
        sf: set(g["row_in_file"].astype(int).tolist())
        for sf, g in meta_test.groupby("source_file")
    }

    for in_path in args.bkg_files:
        stem = Path(in_path).stem
        scored_path = out_dir / f"{stem}_xgb_scored.root"
        allowed = test_rows_by_file.get(in_path, set())
        score_single_root_file(
            model,
            in_path,
            str(scored_path),
            args.tree_name,
            args.branches,
            keep_branches=args.keep_branches,
            allowed_rows=allowed,
        )
        print(f"[xgboost] Wrote scored background file: {scored_path}")

    for in_path in args.sig_files:
        stem = Path(in_path).stem
        scored_path = out_dir / f"{stem}_xgb_scored.root"
        allowed = test_rows_by_file.get(in_path, set())
        score_single_root_file(
            model,
            in_path,
            str(scored_path),
            args.tree_name,
            args.branches,
            keep_branches=args.keep_branches,
            allowed_rows=allowed,
        )
        print(f"[xgboost] Wrote scored signal file: {scored_path}")

    for in_path in data_files:
        stem = Path(in_path).stem
        scored_path = out_dir / f"{stem}_xgb_scored.root"
        score_single_root_file(
            model,
            in_path,
            str(scored_path),
            args.tree_name,
            args.branches,
            keep_branches=args.keep_branches,
        )
        print(f"[xgboost] Wrote scored data file: {scored_path}")

    for in_path in args.det_var_files:
        stem = Path(in_path).stem
        scored_path = out_dir / f"{stem}_xgb_scored.root"
        score_single_root_file(
            model,
            in_path,
            str(scored_path),
            args.tree_name,
            args.branches,
            keep_branches=args.keep_branches,
        )
        print(f"[xgboost] Wrote scored detector-variation file: {scored_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
