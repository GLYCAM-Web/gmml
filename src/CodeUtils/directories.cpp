#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <sys/param.h> // for MIN function
#include <cstring>     // strlen
#include <filesystem>

#include <iostream>

std::string codeUtils::Find_Program_Installation_Directory()
{ // A way to get the program name plus working directory
    char processID[32];
    char pBuffer[256];
    size_t len = sizeof(pBuffer);
    sprintf(processID, "/proc/%d/exe", getpid());
    int bytes = MIN(readlink(processID, pBuffer, len), len - 1);
    if (bytes >= 0)
    {
        pBuffer[bytes] = '\0';
        // std::cout << "processID:" << processID << " pBuffer:" << pBuffer << " bytes:" << bytes << std::endl;
        return codeUtils::SplitFilename(pBuffer);
    }
    return "Error";
}

std::string codeUtils::Find_Program_workingDirectory()
{
    char cCurrentPath[FILENAME_MAX];
    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
    {
        return "Error reading working directory";
    }
    cCurrentPath[strlen(cCurrentPath)] = '/';  // Add a / at the end.
    cCurrentPath[strlen(cCurrentPath)] = '\0'; // Above overwrites the null, the null is important. Respect the nu
    return cCurrentPath;
}

bool codeUtils::doesDirectoryExist(const std::string& pathName)
{
    struct stat info;
    if (stat(pathName.c_str(), &info) != 0)
    {
        return false; // printf( "cannot access %s\n", pathname );
    }
    else if (info.st_mode & S_IFDIR)
    {                // S_ISDIR() doesn't exist on my windows
        return true; // printf( "%s is a directory\n", pathname );
    }
    else
    {
        return false; // printf( "%s is no directory\n", pathname );
    }
}

void codeUtils::ensureDirectoryExists(const std::string& pathName)
{
    if (!codeUtils::doesDirectoryExist(pathName))
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Directory " + pathName + " does not exist");
        throw "Directory " + pathName + " does not exist";
    }
}

std::string codeUtils::getEnvVar(const std::string& key)
{
    char* val = std::getenv(key.c_str());
    return val == NULL ? std::string("") : std::string(val);
}

std::string codeUtils::getGemsHomeDir()
{

    std::filesystem::path directoriesCPPFilePath = __FILE__;
    // TODO: Fix this gross ish
    std::string gemsDirString = directoriesCPPFilePath.parent_path().parent_path().parent_path().parent_path().string();
    gemsDirString             += "/";
    return gemsDirString;
}

std::string codeUtils::getGmmlHomeDir()
{

    std::filesystem::path directoriesCPPFilePath = __FILE__;
    // TODO: Fix this gross ish

    std::string gmmlDirString = directoriesCPPFilePath.parent_path().parent_path().parent_path().string();
    gmmlDirString             += "/";

    std::string gemsDirString = getGemsHomeDir();
    std::string gmmlHome      = gmmlDirString;
    if (gmmlHome.empty())
    {
        std::string gemsHome = gemsDirString;
        if (gemsHome.empty())
        {
            std::string errorMessage =
                "$GMMLHOME and $GEMSHOME environmental variable not set (or std::getenv doesn't work on this system)";
            gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
            throw errorMessage;
        }
        gmmlHome = gemsHome + "/gmml/"; // guessing.
        if (!codeUtils::doesDirectoryExist(gmmlHome))
        {
            std::string errorMessage = "$GMMLHOME not set and directory $GEMSHOME/gmml/ doesn't exist.";
            gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
            throw errorMessage;
        }
    }
    return gmmlHome;
}
