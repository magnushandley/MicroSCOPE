import argparse

import uproot
import pandas as pd
import numpy as np

from typing import List, Optional

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, classification_report

import xgboost as xgb
from pathlib import Path

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
    test_size: float = 0.5,
    random_state: int = 1337,
    balance_sig_to_bkg: bool = True,
):
    """Prepare X/y/weights and split into train/test.

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

    # Build arrays
    X = df[feature_columns]
    y = df["label"].astype(int)
    w = df[weight_column].astype(float)

    # Split (stratify to preserve class fractions)
    X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
        X,
        y,
        w,
        test_size=test_size,
        random_state=random_state,
        stratify=y,
    )

    # Keep metadata to allow filtering back to original files
    meta_cols = ["source_file", "row_in_file", "label", weight_column]
    meta = df[meta_cols]
    meta_train = meta.loc[X_train.index].copy()
    meta_test = meta.loc[X_test.index].copy()

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


# --- Helper functions for scoring and writing ROOT ---
def add_bdt_score(
    booster: xgb.Booster,
    df: pd.DataFrame,
    feature_columns: List[str],
    score_column: str = "bdt_score",
) -> pd.DataFrame:
    """Return a copy of df with an added BDT score column."""
    df_out = df.copy()
    dmat = xgb.DMatrix(df_out[feature_columns])
    df_out[score_column] = booster.predict(dmat)
    return df_out


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
        Mapping of tree_name -> DataFrame to write.
    """

    with uproot.recreate(out_path) as f:
        for tree_name, df in trees.items():
            # Drop non-numeric columns (e.g. source_file paths)
            cols_numeric = [
                c for c in df.columns
                if pd.api.types.is_numeric_dtype(df[c])
            ]
            branch_arrays = {c: df[c].to_numpy() for c in cols_numeric}

            tree_path = [part for part in tree_name.split("/") if part]
            tree_dir = f
            for directory in tree_path[:-1]:
                tree_dir = tree_dir.mkdir(directory)

            tree_leaf_name = tree_path[-1]
            branch_types = {
                branch_name: branch_array.dtype
                for branch_name, branch_array in branch_arrays.items()
            }
            tree_dir.mktree(tree_leaf_name, branch_types)
            tree_dir[tree_leaf_name].extend(branch_arrays)


def score_single_root_file(
    booster: xgb.Booster,
    in_path: str,
    out_path: str,
    tree_name: str,
    feature_columns: List[str],
    allowed_rows: Optional[set] = None,
):
    """Load one ROOT file, add bdt_score branch, and write to a new ROOT file.

    The output file will contain a TTree with the same name as `tree_name`.
    Non-numeric helper columns (e.g. source_file) are dropped automatically.
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

    df_scored = add_bdt_score(booster, df, feature_columns=feature_columns, score_column="bdt_score")

    # Write a single tree with the original tree_name
    write_trees_to_root(out_path, {tree_name: df_scored})

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Train an XGBoost BDT and write scored ROOT files."
    )
    parser.add_argument("--config", help="Optional MicroSCOPE config path passed through by PythonModule.")
    parser.add_argument("--bkg-files", nargs="+", required=True, help="Background ROOT input files.")
    parser.add_argument("--sig-files", nargs="+", required=True, help="Signal ROOT input files.")
    parser.add_argument("--data-file", required=True, help="Data ROOT input file to score.")
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
    parser.add_argument("--output-dir", default=".", help="Directory for model and scored ROOT outputs.")
    parser.add_argument("--model-output", default="xgb_bdt.json", help="Model filename or path.")
    parser.add_argument("--test-size", type=float, default=0.5, help="Held-out test fraction.")
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
    if not args.data_file:
        raise ValueError("--data-file cannot be empty")
    if not args.tree_name:
        raise ValueError("--tree-name cannot be empty")
    if not args.branches:
        raise ValueError("--branches must contain at least one branch")
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
            allowed_rows=allowed,
        )
        print(f"[xgboost] Wrote scored signal file: {scored_path}")

    stem = Path(args.data_file).stem
    scored_path = out_dir / f"{stem}_xgb_scored.root"
    score_single_root_file(model, args.data_file, str(scored_path), args.tree_name, args.branches)
    print(f"[xgboost] Wrote scored data file: {scored_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
