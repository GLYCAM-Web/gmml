#include <cstdlib> // getenv
#include <sys/stat.h> // To check if file exists using stat
#include <sys/types.h> // The s_IFDIR
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"

bool codeutils::doesFileExist(const std::string& fileName)
{
    struct stat buffer;
    return (stat (fileName.c_str(), &buffer) == 0);
}

void codeutils::ensureFileExists(const std::string& fileName)
{
    if (!codeutils::doesFileExist(fileName))
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "File " + fileName + " does not exist");
        throw "File " + fileName + " does not exist";
    }
}

bool codeutils::doesDirectoryExist(const std::string& pathName)
{
    struct stat info;
    if( stat( pathName.c_str(), &info ) != 0 )
        return false; //printf( "cannot access %s\n", pathname );
    else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
        return true; //printf( "%s is a directory\n", pathname );
    else
        return false; //printf( "%s is no directory\n", pathname );
}

void codeutils::ensureDirectoryExists (const std::string& pathName)
{
    if (!codeutils::doesDirectoryExist(pathName))
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Directory " + pathName + " does not exist");
        throw "Directory " + pathName + " does not exist";
    }
}

std::string codeutils::getEnvVar( std::string const & key )
{
    char * val = std::getenv( key.c_str() );
    return val == NULL ? std::string("") : std::string(val);
}

std::string codeutils::getGmmlHomeDir()
{
    std::string gmmlHome = getEnvVar("GMMLHOME");
    if (gmmlHome.empty())
    {
        std::string gemsHome = getEnvVar("GEMSHOME");
        if (gemsHome.empty())
        {
            std::string errorMessage = "$GMMLHOME and $GEMSHOME environmental variable not set (or std::getenv doesn't work on this system)";
            gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
            throw errorMessage;
        }
        gmmlHome = gemsHome + "/gmml/"; // guessing.
        if(!codeutils::doesDirectoryExist(gmmlHome))
        {
            std::string errorMessage = "$GMMLHOME not set and directory $GEMSHOME/gmml/ doesn't exist.";
            gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
            throw errorMessage;
        }
    }
    return gmmlHome;
}
