#include "Modules/PythonModule.hxx"

#include "Utils/ConfigUtils.hxx"

#include <cerrno>
#include <cctype>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>
#include <utility>

using namespace Analysis;

PythonModule::PythonModule(const TEnv& cfg, std::string configPath)
    : Module(cfg)
    , fConfigPath(std::move(configPath))
{
    ParseConfig();
}

void PythonModule::ParseConfig()
{
    fScript = RequireValue(Cfg(), "Python.Script");
    fPythonExecutable = Cfg().GetValue("Python.PythonExecutable", "python3");
    fWorkingDirectory = Cfg().GetValue("Python.WorkingDirectory", ".");
    fPassConfig = Cfg().GetValue("Python.PassConfig", false);
    fConfigArg = Cfg().GetValue("Python.ConfigArg", "--config");

    if (fScript.empty()) {
        throw std::runtime_error("[PythonModule] Python.Script cannot be empty");
    }
    if (fPythonExecutable.empty()) {
        throw std::runtime_error("[PythonModule] Python.PythonExecutable cannot be empty");
    }
    if (fPassConfig && fConfigArg.empty()) {
        throw std::runtime_error("[PythonModule] Python.ConfigArg cannot be empty when PassConfig is true");
    }

    fCommand = BuildCommand();
}

std::vector<std::string> PythonModule::BuildCommand() const
{
    std::vector<std::string> command;
    command.push_back(fPythonExecutable);
    command.push_back(fScript);

    if (fPassConfig) {
        command.push_back(fConfigArg);
        command.push_back(fConfigPath);
    }

    const auto argumentNames = SplitWhitespace(Cfg().GetValue("Python.Arguments", ""));

    for (const auto& argumentName : argumentNames) {
        const std::string valueKey = "Python.Arg." + argumentName;
        const std::string value = RequireValue(Cfg(), valueKey);
        command.push_back(ToCliFlag(argumentName));

        const auto tokens = SplitWhitespace(value);
        command.insert(command.end(), tokens.begin(), tokens.end());
    }

    return command;
}

void PythonModule::Initialise()
{
    std::cout << "[PythonModule] Executable: " << fPythonExecutable << "\n";
    std::cout << "[PythonModule] Script: " << fScript << "\n";
    std::cout << "[PythonModule] WorkingDirectory: " << fWorkingDirectory << "\n";
    std::cout << "[PythonModule] Command: " << FormatCommand(fCommand) << "\n";

    const int exitCode = RunCommand();
    if (exitCode != 0) {
        std::ostringstream message;
        message << "[PythonModule] Python process failed with exit code " << exitCode;
        throw std::runtime_error(message.str());
    }
}

int PythonModule::RunCommand() const
{
    int pipefd[2];
    if (pipe(pipefd) == -1) {
        throw std::runtime_error(std::string("[PythonModule] pipe failed: ")
                                 + std::strerror(errno));
    }

    const pid_t pid = fork();
    if (pid == -1) {
        close(pipefd[0]);
        close(pipefd[1]);
        throw std::runtime_error(std::string("[PythonModule] fork failed: ")
                                 + std::strerror(errno));
    }

    if (pid == 0) {
        close(pipefd[0]);
        dup2(pipefd[1], STDOUT_FILENO);
        dup2(pipefd[1], STDERR_FILENO);
        close(pipefd[1]);

        if (!fWorkingDirectory.empty() && chdir(fWorkingDirectory.c_str()) != 0) {
            std::cerr << "chdir(" << fWorkingDirectory << ") failed: "
                      << std::strerror(errno) << "\n";
            _exit(127);
        }

        std::vector<char*> argv;
        argv.reserve(fCommand.size() + 1);
        for (const auto& part : fCommand) {
            argv.push_back(const_cast<char*>(part.c_str()));
        }
        argv.push_back(nullptr);

        execvp(argv[0], argv.data());
        std::cerr << "execvp(" << argv[0] << ") failed: " << std::strerror(errno) << "\n";
        _exit(127);
    }

    close(pipefd[1]);
    StreamWithPrefix(pipefd[0], "[PythonModule] ");
    close(pipefd[0]);

    int status = 0;
    while (waitpid(pid, &status, 0) == -1) {
        if (errno != EINTR) {
            throw std::runtime_error(std::string("[PythonModule] waitpid failed: ")
                                     + std::strerror(errno));
        }
    }

    if (WIFEXITED(status)) {
        return WEXITSTATUS(status);
    }

    if (WIFSIGNALED(status)) {
        std::ostringstream message;
        message << "[PythonModule] Python process terminated by signal "
                << WTERMSIG(status);
        throw std::runtime_error(message.str());
    }

    throw std::runtime_error("[PythonModule] Python process ended unexpectedly");
}

std::vector<std::string> PythonModule::SplitWhitespace(const std::string& value)
{
    std::vector<std::string> tokens;
    std::stringstream ss(value);
    std::string token;
    while (ss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string PythonModule::RequireValue(const TEnv& cfg, const std::string& key)
{
    if (!ConfigHasKey(cfg, key)) {
        throw std::runtime_error("[PythonModule] Missing required config key: " + key);
    }
    return cfg.GetValue(key.c_str(), "");
}

std::string PythonModule::ToCliFlag(const std::string& name)
{
    std::string flag = "--";
    bool lastWasDash = false;

    for (std::size_t i = 0; i < name.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(name[i]);
        if (name[i] == '_' || name[i] == '-' || name[i] == '.' || std::isspace(c)) {
            if (!lastWasDash && flag.size() > 2) {
                flag.push_back('-');
                lastWasDash = true;
            }
            continue;
        }

        const bool isUpper = std::isupper(c) != 0;
        const bool previousIsLowerOrDigit =
            i > 0 && (std::islower(static_cast<unsigned char>(name[i - 1])) != 0
                      || std::isdigit(static_cast<unsigned char>(name[i - 1])) != 0);
        const bool nextIsLower =
            i + 1 < name.size()
            && std::islower(static_cast<unsigned char>(name[i + 1])) != 0;
        const bool previousIsUpper =
            i > 0 && std::isupper(static_cast<unsigned char>(name[i - 1])) != 0;

        if (isUpper && flag.size() > 2
            && (previousIsLowerOrDigit || (previousIsUpper && nextIsLower))
            && !lastWasDash) {
            flag.push_back('-');
        }

        flag.push_back(static_cast<char>(std::tolower(c)));
        lastWasDash = false;
    }

    if (!flag.empty() && flag.back() == '-') {
        flag.pop_back();
    }
    return flag;
}

std::string PythonModule::FormatCommand(const std::vector<std::string>& command)
{
    std::ostringstream out;
    for (std::size_t i = 0; i < command.size(); ++i) {
        if (i > 0) out << ' ';

        const bool needsQuotes = command[i].find_first_of(" \t\n'\"\\") != std::string::npos;
        if (!needsQuotes) {
            out << command[i];
            continue;
        }

        out << '\'';
        for (const char c : command[i]) {
            if (c == '\'') {
                out << "'\\''";
            } else {
                out << c;
            }
        }
        out << '\'';
    }
    return out.str();
}

void PythonModule::StreamWithPrefix(int fd, const std::string& prefix)
{
    std::string pending;
    char buffer[4096];

    while (true) {
        const ssize_t nRead = read(fd, buffer, sizeof(buffer));
        if (nRead == 0) break;
        if (nRead < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error("[PythonModule] read failed: " + std::string(std::strerror(errno)));
        }

        pending.append(buffer, static_cast<std::size_t>(nRead));
        std::size_t lineStart = 0;
        while (true) {
            const std::size_t newline = pending.find('\n', lineStart);
            if (newline == std::string::npos) break;

            std::cout << prefix << pending.substr(lineStart, newline - lineStart) << "\n";
            lineStart = newline + 1;
        }
        pending.erase(0, lineStart);
    }

    if (!pending.empty()) {
        std::cout << prefix << pending << "\n";
    }
}
