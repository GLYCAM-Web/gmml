#ifndef GMML_INCLUDES_CODEUTILS_FILES_HPP
#define GMML_INCLUDES_CODEUTILS_FILES_HPP
#include <sys/stat.h> // To check if file exists using stat
#include <sys/types.h> // The s_IFDIR
#include <string>
#include "includes/CodeUtils/logging.hpp"

namespace gmml
{
inline bool doesFileExist(const std::string& fileName)
{
    struct stat buffer;
    return (stat (fileName.c_str(), &buffer) == 0);
}

inline void ensureFileExists(const std::string& fileName)
{
    if (!doesFileExist(fileName))
    {

        throw "File " + fileName + " does not exist";
    }
}

inline bool doesDirectoryExist(const std::string& pathName)
{
    struct stat info;
    if( stat( pathName.c_str(), &info ) != 0 )
        return false; //printf( "cannot access %s\n", pathname );
    else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
        return true; //printf( "%s is a directory\n", pathname );
    else
        return false; //printf( "%s is no directory\n", pathname );
}

inline void ensureDirectoryExists (const std::string& pathName)
{
    if (!doesDirectoryExist(pathName))
    {
        throw "Directory " + pathName + " does not exist";
    }
}


}
#endif // GMML_INCLUDES_CODEUTILS_FILES_HPP
