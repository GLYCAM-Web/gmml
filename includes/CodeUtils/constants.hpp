#ifndef GMML_INCLUDES_CODEUTILS_CONSTANTS_HPP
#define GMML_INCLUDES_CODEUTILS_CONSTANTS_HPP
#include <math.h> // atan

namespace constants
{

const double DEFAULT_ANGLE = 109.4;
const double CARBON_SURFACE_AREA = 36.31681103;
const double MAX_RESIDUE_DIAMETER = 9.0;
const double PI_RADIAN = 4.0*atan(1.0);
const double PI_DEGREE = 180.0;
const double dNotSet = 123456789.0;
const int iNotSet = -123456;
const std::string sNotSet = "?";

}// namespace

#endif
