#ifndef LAZY_PRINTERS_HPP
#define LAZY_PRINTERS_HPP

#include <iostream>

// LAZY INFORMATION OUTPUT
inline void lazyInfo(int lineCalled_t, const char *funcName_t)
{
  std::cout << "****************************************\n" << "\tINFORMATION: " << funcName_t << "\n" << "****************************************\n";
  std::cout << "Function Name: " << funcName_t << std::endl;
  std::cout << "Line Number: " << lineCalled_t << std::endl << std::endl;
}

inline void lazyInfo(int lineCalled_t, const char *funcName_t, std::string infoToPass_t)
{
  std::cout << "****************************************\n" << "\t\tINFORMATION\n" << "****************************************\n";
  std::cout << "Function Name: " << funcName_t << std::endl;
  std::cout << "Line Number: " << lineCalled_t << std::endl;
  std::cout << "Msg: " << infoToPass_t << std::endl << std::endl;
}

/* NOTE: Please note that these are to temporarily keep track of unwanted behavior
 * 			and we will eventually replace these with correct error handling etc.
 * 			but while the program and logic are in flux I would prefer to keep
 * 			our error handling ambiguous.
 *
 */
inline void badBehavior(int lineBroke_t, const char *funcNameBroke_t)
{
  std::cout << "****************************************\n" << "\t\tBORKED\n" << "****************************************\n";
  std::cout << "Borked Function Name: " << funcNameBroke_t << std::endl;
  std::cout << "Line Number: " << lineBroke_t << std::endl << std::endl;
}

inline void badBehavior(int lineBroke_t, const char *funcNameBroke_t, std::string infoToPass_t)
{
  std::cout << "****************************************\n" << "\t\tBORKED\n" << "****************************************\n";
  std::cout << "Borked Function Name: " << funcNameBroke_t << std::endl;
  std::cout << "Line Number: " << lineBroke_t << std::endl << std::endl;
  std::cout << "Msg: " << infoToPass_t << std::endl << std::endl;
}

#endif // LAZY_PRINTERS_HPP
