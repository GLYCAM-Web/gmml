#ifndef PARAMETERFILEDIHEDRALTERM_HPP
#define PARAMETERFILEDIHEDRALTERM_HPP

#include <string>
#include <iostream>

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

///Delaram
            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to access to the factor attribute of the current object
              * The attribute is set by the contents of the given file
              * @return factor_ of the current object of this class
              */
            double GetFactor();
            /*! \fn
              * An accessor function in order to access to force constant attribute of the current object
              * The attribute is set by the contents of the given file
              * @return force_constant_ of the current object of this class
              */
            double GetForceConstant();
            /*! \fn
              * An accessor function in order to access to phase attribute of the current object
              * The attribute is set by the contents of the given file
              * @return phase_ of the current object of this class
              */
            double GetPhase();
            /*! \fn
              * An accessor function in order to access to periodicity attribute of the current object
              * The attribute is set by the contents of the given file
              * @return periodicity_ of the current object of this class
              */
            double GetPeriodicity();
            /*! \fn
              * An accessor function in order to access to dscr attribute of the current object
              * @return dscr_ attribute of the current object of this class
              */
            std::string GetDscr();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the factor of the current object
              * Set the factor_ attribute of the current term
              * @param factor The factor of the current object
              */
            void SetFactor(double factor);
            /*! \fn
              * A mutator function in order to set the force_constant of the current object
              * Set the force_constant_ attribute of the current term
              * @param force_constant The force_constant of the current object
              */
            void SetForceConstant(double force_constant);
            /*! \fn
              * A mutator function in order to set the phase of the current object
              * Set the phase_ attribute of the current term
              * @param phase The phase of the current object
              */
            void SetPhase(double phase);
            /*! \fn
              * A mutator function in order to set the periodicity of the current object
              * Set the periodicity_ attribute of the current term
              * @param periodicity The periodicity of the current object
              */
            void SetPeriodicity(double periodicity);
            /*! \fn
              * A mutator function in order to set the dscr of the current object
              * Set the dscr_ attribute of the current term
              * @param dscr The dscr attribute of the current object
              */
            void SetDscr( std::string dscr);
///Delaram

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
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
