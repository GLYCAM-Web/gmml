#ifndef PARAMETERFILEDIHEDRALTERM_HPP
#define PARAMETERFILEDIHEDRALTERM_HPP

#include <string>

namespace ParameterFileSpace
{
    class ParameterFileDihedralTerm
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            ParameterFileDihedralTerm();
            /*! \fn
              * Constructor with essential parameters
              * @param factor
              * @param force_constant
              * @param phase
              * @param periodicity
              * @param dscr
              */
            ParameterFileDihedralTerm(double factor, double force_constant, double phase, double periodicity, const std::string& dscr = "");

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
            double factor_;                         /*!< Division factor; Fill by the second column of the 5th section of the parameter file (6th section doesn't have this value)*/
            double force_constant_;                 /*!< Dihedral force constant; Fill by the third column of the 5th section or second column of the 6th section of the parameter file*/
            double phase_;                          /*!< Dihedral phase; Fill by the fourth column of the 5th section of third column of the 6th section of the parameter file*/
            double periodicity_;                    /*!< Dihedral periodicity; Fill by the 5th column of the 5th section or fourth column of the 6th section of the parameter file*/
            std::string dscr_;                      /*!< Dihedral description; Fill by the last column of the 5th or 6th section of the parameter file*/
    };
}

#endif // PARAMETERFILEDIHEDRALTERM_HPP
