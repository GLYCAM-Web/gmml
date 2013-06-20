#ifndef PARAMETERFILEANGLE_HPP
#define PARAMETERFILEANGLE_HPP

#include <string>
#include <vector>

namespace ParameterFileSpace
{
    class ParameterFileAngle
    {
        public:
            /////////////////////////////// CONSTRUCTOR ////////////////////////////////
            ParameterFileAngle();
            ParameterFileAngle(const std::vector<std::string>& types, double force_constant, double angle, const std::string& dscr = "");

            /////////////////////////// DISPLAY FUNCTION ///////////////////////////////
            void Print(std::ostream& out);

            /////////////////////////////// ATTRIBUTES /////////////////////////////////
            std::vector<std::string> types_;        // Atom types involved in an angle; Fill by the firs column of the fourth section of the parameter file split by '-'
            double force_constant_;                 // Angle force constant; Fill by the second column of the fourth section of the parameter file
            double angle_;                          // Angle value; Fill by the third column third column of the fourth section of the parameter file
            std::string dscr_;                      // Angle description; Fill by the fourth column of the fourth section of the parameter file
            // A sample of the fourth (Angle) sectio of the parameter file: C -C -O     80.0      120.00    Junmei et al, 1999 acrolein
    };
}

#endif // PARAMETERFILEANGLE_HPP
