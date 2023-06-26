#include <cstdlib>     // getenv
#include <sys/stat.h>  // To check if file exists using stat
#include <sys/types.h> // The s_IFDIR
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"

bool codeUtils::doesFileExist(const std::string& fileName)
{
    struct stat buffer;
    return (stat(fileName.c_str(), &buffer) == 0);
}

void codeUtils::ensureFileExists(const std::string& fileName)
{
    if (!codeUtils::doesFileExist(fileName))
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "File " + fileName + " does not exist");
        throw "File " + fileName + " does not exist";
    }
}

std::string codeUtils::SplitFilename(const std::string& str)
{
    // std::cout << "Splitting: " << str << '\n';
    std::size_t found = str.find_last_of("/\\");
    //  std::cout << " path: " << str.substr(0,found) << '\n';
    //  std::cout << " file: " << str.substr(found+1) << '\n';
    return str.substr(0, found);
}
