#ifndef GMML_INCLUDES_CODEUTILS_FILES_HPP
#define GMML_INCLUDES_CODEUTILS_FILES_HPP

#include <string>

namespace codeutils
{
bool doesFileExist(const std::string& fileName);
void ensureFileExists(const std::string& fileName);
bool doesDirectoryExist(const std::string& pathName);
void ensureDirectoryExists (const std::string& pathName);
std::string getEnvVar( std::string const & key );
std::string getGmmlHomeDir();
}
#endif //GMML_INCLUDES_CODEUTILS_FILES_HPP
