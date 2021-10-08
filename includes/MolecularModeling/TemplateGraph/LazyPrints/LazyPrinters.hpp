#ifndef LAZY_PRINTERS_HPP
#define LAZY_PRINTERS_HPP

#include <iostream>

// LAZY INFORMATION OUTPUT
inline void lazyInfo(int lineCalled_t, const char *funcName_t)
{
  std::stringstream logss;
  logss << "****************************************\n" << "\tINFORMATION: " << funcName_t << "\n" << "****************************************\n";
  logss << "Function Name: " << funcName_t << std::endl;
  logss << "Line Number: " << lineCalled_t << std::endl << std::endl;
  gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
}

inline void lazyInfo(int lineCalled_t, const char *funcName_t, std::string infoToPass_t)
{
  std::stringstream logss;
  logss << "****************************************\n" << "\t\tINFORMATION\n" << "****************************************\n";
  logss << "Function Name: " << funcName_t << std::endl;
  logss << "Line Number: " << lineCalled_t << std::endl;
  logss << "Msg: " << infoToPass_t << std::endl << std::endl;
  gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
}

/* NOTE: Please note that these are to temporarily keep track of unwanted behavior
 * 			and we will eventually replace these with correct error handling etc.
 * 			but while the program and logic are in flux I would prefer to keep
 * 			our error handling ambiguous.
 *
 */
inline void badBehavior(int lineBroke_t, const char *funcNameBroke_t)
{
  std::stringstream logss;
  logss << "****************************************\n" << "\t\tBORKED\n" << "****************************************\n";
  logss << "Borked Function Name: " << funcNameBroke_t << std::endl;
  logss << "Line Number: " << lineBroke_t << std::endl << std::endl;
  gmml::log(__LINE__, __FILE__, gmml::ERR, logss.str());
}

inline void badBehavior(int lineBroke_t, const char *funcNameBroke_t, std::string infoToPass_t)
{
  std::stringstream logss;
  logss << "****************************************\n" << "\t\tBORKED\n" << "****************************************\n";
  logss << "Borked Function Name: " << funcNameBroke_t << std::endl;
  logss << "Line Number: " << lineBroke_t << std::endl << std::endl;
  logss << "Msg: " << infoToPass_t << std::endl << std::endl;
  gmml::log(__LINE__, __FILE__, gmml::ERR, logss.str());
}

#endif // LAZY_PRINTERS_HPP
