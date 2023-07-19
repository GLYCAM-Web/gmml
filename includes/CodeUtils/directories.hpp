#ifndef GMML_INCLUDES_CODEUTILS_DIRECTORIES_HPP
#define GMML_INCLUDES_CODEUTILS_DIRECTORIES_HPP
#include <sys/stat.h>  // To check if file exists using stat
#include <sys/types.h> // The s_IFDIR
#include <string>

// #include "includes/CodeUtils/logging.hpp"

namespace codeUtils
{

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

    std::string Find_Program_Installation_Directory();
    std::string Find_Program_workingDirectory();
    bool doesDirectoryExist(const std::string& pathName);
    void ensureDirectoryExists(const std::string& pathName);
    std::string getEnvVar(const std::string& key);
    std::string getGmmlHomeDir();
    std::string getGemsHomeDir();
} // namespace codeUtils
#endif // GMML_INCLUDES_CODEUTILS_DIRECTORIES_HPP
