#ifndef GMML_INCLUDES_CODEUTILS_LOGGING_HPP
#define GMML_INCLUDES_CODEUTILS_LOGGING_HPP

#include <string>

namespace gmml
{
  /* \todo DM 1-21-22
  add __func__ to show the function name in the log
  Add or edit logLevel values
  */
    enum LogLevel
    {
        INF,
        ERR,
        WAR
    };
/*! \brief A function to write information/warning/error messages into a log file

      * @pre GMMLLOG environment variable has to be set, "GMML_Log.txt" will be ignored by git
      * @param line The line number producing the message
      * @param file_path The file path of the file which the message has been produced within in
      * @param level The type of the produced message INF/WAR/ERR
      * @param msg The message content that has been produced
      * @param out_file_name The name of the output log file
      */
    void log(int line, std::string file_path, LogLevel level, std::string msg, std::string out_file_name = "");
}
#endif //GMML_INCLUDES_CODEUTILS_LOGGING_HPP
