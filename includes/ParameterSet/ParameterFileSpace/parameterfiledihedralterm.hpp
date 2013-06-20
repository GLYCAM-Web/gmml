#ifndef PARAMETERFILEDIHEDRALTERM_HPP
#define PARAMETERFILEDIHEDRALTERM_HPP

#include <string>

namespace ParameterFileSpace
{
    class ParameterFileDihedralTerm
    {
        public:
            /////////////////////////////// CONSTRUCTOR ////////////////////////////////
            ParameterFileDihedralTerm();
            ParameterFileDihedralTerm(double factor, double force_constant, double phase, double periodicity, const std::string& dscr = "");

            //////////////////////////// DISPLAY FUNCTION //////////////////////////////
            void Print(std::ostream& out);

            /////////////////////////////// ATTRIBUTES /////////////////////////////////
            double factor_;                         // Division factor; Fill by the second column of the 5th section of the parameter file (6th section doesn't have this value)
            double force_constant_;                 // Dihedral force constant; Fill by the third column of the 5th section or second column of the 6th section of the parameter file
            double phase_;                          // Dihedral phase; Fill by the fourth column of the 5th section of third column of the 6th section of the parameter file
            double periodicity_;                    // Dihedral periodicity; Fill by the 5th column of the 5th section or fourth column of the 6th section of the parameter file
            std::string dscr_;                      // Dihedral description; Fill by the last column of the 5th or 6th section of the parameter file
    };
}

#endif // PARAMETERFILEDIHEDRALTERM_HPP
