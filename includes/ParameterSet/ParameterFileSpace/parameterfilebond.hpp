#ifndef PARAMETERFILEBOND_HPP
#define PARAMETERFILEBOND_HPP

#include <vector>
#include <string>

namespace ParameterFileSpace
{
    class ParameterFileBond
    {
        public:
            /////////////////////////////// CONSTRUCTOR ////////////////////////////////
            ParameterFileBond();
            ParameterFileBond(const std::vector<std::string>& types, double force_constant, double length, const std::string& dscr = "");
            ParameterFileBond(const std::vector<std::string>& types, double force_constant, double length, const std::vector<double>& hbond_coefficients,
                              const std::string& dscr = "");

            /////////////////////////// DISPLAY FUNCTION //////////////////////////////
            void Print(std::ostream& out);

            ///////////////////////////// ATTRIBUTES //////////////////////////////////
            std::vector<std::string> types_;        // Atom types involved in a bond; Fill by the first column of the third section of the parameter file split by '-'
            double force_constant_;                 // Bond force constant; Fill by the second column of the third section of the parameter file
            double length_;                         // Bond Length; Fill by the third column of the third section of the parameter file
            std::string dscr_;                      // Bond description; Fill by the fourth column of the third section of the parameter file
            // A sample of the third (Bond) section of the parameter file: C -C   310.0    1.525       Junmei et al, 1999

            std::vector<double> hbond_coefficients_;  // Hydrogen bond coefficients; Fill by the 7th section of the parameter file
            // A sample of the 7th (hydrogen bond) section of the parameter file:   HW  OW  0000.     0000.                                4.  flag for fast water
    };
}

#endif // PARAMETERFILEBOND_HPP
