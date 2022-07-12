#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/files.hpp"


#include <sys/param.h> // for MIN function
#include <cstring> // strlen


std::string gmml::Find_Program_Installation_Directory()
{   // A way to get the program name plus working directory
    char processID[32];
    char pBuffer[256];
    size_t len = sizeof(pBuffer);
    sprintf(processID, "/proc/%d/exe", getpid());
    int bytes = MIN(readlink(processID, pBuffer, len), len - 1);
    if(bytes >= 0)
    {
        pBuffer[bytes] = '\0';
        //std::cout << "processID:" << processID << " pBuffer:" << pBuffer << " bytes:" << bytes << std::endl;
        return gmml::SplitFilename (pBuffer);
    }
    return "Error";
}

std::string gmml::Find_Program_workingDirectory()
{
    char cCurrentPath[FILENAME_MAX];
    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
        return "Error reading working directory";
    cCurrentPath[strlen(cCurrentPath)] = '/'; // Add a / at the end.
    cCurrentPath[strlen(cCurrentPath)] = '\0'; // Above overwrites the null, the null is important. Respect the nu
    return cCurrentPath;
}
