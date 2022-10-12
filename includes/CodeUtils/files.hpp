#ifndef GMML_INCLUDES_CODEUTILS_FILES_HPP
#define GMML_INCLUDES_CODEUTILS_FILES_HPP
#include <sys/stat.h> // To check if file exists using stat
#include <sys/types.h> // The s_IFDIR
#include <string>
//#include "includes/CodeUtils/logging.hpp"

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


inline std::string SplitFilename (const std::string& str)
{
	//std::cout << "Splitting: " << str << '\n';
	std::size_t found = str.find_last_of("/\\");
	//  std::cout << " path: " << str.substr(0,found) << '\n';
	//  std::cout << " file: " << str.substr(found+1) << '\n';
	return str.substr(0,found);
}


}
#endif //GMML_INCLUDES_CODEUTILS_FILES_HPP
