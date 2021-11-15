#ifndef GMML_INCLUDES_INTERNALPROGRAMS_IO_HPP
#define GMML_INCLUDES_INTERNALPROGRAMS_IO_HPP
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <fstream>      // std::ifstream
#include <stdio.h>  /* defines FILENAME_MAX */
#include <sys/param.h> // for MIN function

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

bool fileExists(std::string fileName);
std::string SplitFilename (const std::string& str);
std::string Find_Program_Installation_Directory();
std::string Find_Program_workingDirectory();

template<typename Out>
inline void split(const std::string &s, char delim, Out result)
{
	std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) 
    {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim);
#endif // GMML_INCLUDES_INTERNALPROGRAMS_IO_HPP
