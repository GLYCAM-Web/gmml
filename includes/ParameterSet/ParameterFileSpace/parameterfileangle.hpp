#ifndef PARAMETERFILEANGLE_HPP
#define PARAMETERFILEANGLE_HPP

#include <string>
#include <vector>

namespace ParameterFileSpace
{
    class ParameterFileAngle
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            ParameterFileAngle();
            /*! \fn
              * Constructor with required parameters
              * @param types A vector of three atom types involving in an angle
              * @param force_constant A double value of constant force between three atoms in an angle
              * @param angle The actual value of angle between three atoms
              * @param dscr A short description for each angle mentioned in a parameter file
              */
            ParameterFileAngle(const std::vector<std::string>& types, double force_constant, double angle, const std::string& dscr = "");

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::vector<std::string> types_;        /*!< Atom types involved in an angle; Fill by the firs column of the fourth section of the parameter file split by '-'*/
            double force_constant_;                 /*!< Angle force constant; Fill by the second column of the fourth section of the parameter file*/
            double angle_;                          /*!< Angle value; Fill by the third column third column of the fourth section of the parameter file*/
            std::string dscr_;                      /*!< Angle description; Fill by the fourth column of the fourth section of the parameter file*/
            /*!< A sample of the fourth (Angle) sections of the parameter file: C -C -O     80.0      120.00    Junmei et al, 1999 acrolein*/
    };
}

#endif // PARAMETERFILEANGLE_HPP
