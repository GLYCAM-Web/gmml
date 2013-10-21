#ifndef PARAMETERFILEBOND_HPP
#define PARAMETERFILEBOND_HPP

#include <vector>
#include <string>

namespace ParameterFileSpace
{
    class ParameterFileBond
    {
        public:
            ///////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            ParameterFileBond();
            /*! \fn
              * Constructor with essential parameters
              * @param types Vector of strings of atom types that are involved in a bond
              * @param force_constant Value of constant force that is established between atoms involving in a bond
              * @param length Length between two atoms involving in a bond
              * @param dscr A description that is mentioned in the bond section of the parameter file
              */
            ParameterFileBond(const std::vector<std::string>& types, double force_constant, double length, const std::string& dscr = "");
            /*! \fn
              * Constructor with essential parameters
              * @param types Vector of strings of atom types that are involved in a bond
              * @param force_constant Value of constant force that is established between atoms involving in a bond
              * @param length Length between two atoms involving in a bond
              * @param hbond_coefficient Two essential coefficients for different hydrogen bonds
              * @param dscr A description that is mentioned in the bond section of the parameter file
              */
            ParameterFileBond(const std::vector<std::string>& types, double force_constant, double length, const std::vector<double>& hbond_coefficients,
                              const std::string& dscr = "");

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out);

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::vector<std::string> types_;        /*!< Atom types involved in a bond; Fill by the first column of the third section of the parameter file split by '-'*/
            double force_constant_;                 /*!< Bond force constant; Fill by the second column of the third section of the parameter file*/
            double length_;                         /*!< Bond Length; Fill by the third column of the third section of the parameter file*/
            std::string dscr_;                      /*!< Bond description; Fill by the fourth column of the third section of the parameter file*/
            /*!< A sample of the third (Bond) section of the parameter file: C -C   310.0    1.525       Junmei et al, 1999*/

            std::vector<double> hbond_coefficients_;  /*!< Hydrogen bond coefficients; Fill by the 7th section of the parameter file*/
            /*!< A sample of the 7th (hydrogen bond) section of the parameter file:   HW  OW  0000.     0000.                                4.  flag for fast water*/
    };
}

#endif // PARAMETERFILEBOND_HPP
