#ifndef GMML_INCLUDES_CODEUTILS_LOGGING_HPP
#define GMML_INCLUDES_CODEUTILS_LOGGING_HPP

#include <string>

namespace gmml
{
    enum LogLevel
    {
        INF,
        ERR,
        WAR
    };
/*! \fn
      * A function in order to write the information/warning/error messages produced by the program into a log file
      * @param line The line number producing the message
      * @param file_path The file path of the file which the message has been produced within in
      * @param level The type of the produced message INF/WAR/ERR
      * @param msg The message content that has been produced
      * @param out_file_name The name of the output log file
      */
    void log(int line, std::string file_path, LogLevel level, std::string msg, std::string out_file_name = "");
}
#endif //GMML_INCLUDES_CODEUTILS_LOGGING_HPP
