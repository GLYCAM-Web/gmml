#ifndef PARAMETERFILEBOND_HPP
#define PARAMETERFILEBOND_HPP

#include <vector>
#include <string>
#include <iostream>

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
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the atom types involved in an angle in the current object
              * @return types_ attribute of the current object of this class
              */
            std::vector<std::string> GetTypes();
            /*! \fn
              * An accessor function in order to access to access to the force constant attribute of the current object
              * The attribute is set by the contents of the given file
              * @return force_constant_ of the current object of this class
              */
            double GetForceConstant();
            /*! \fn
              * An accessor function in order to access to length attribute of the current object
              * The attribute is set by the contents of the given file
              * @return length_ of the current object of this class
              */
            double GetLength();
            /*! \fn
              * An accessor function in order to access to angle descripion of the current object
              * @return dscr_ attribute of the current object of this class
              */
            std::string GetDscr();
            /*! \fn
              * An accessor function in order to access to the hbond coefficients involved in an angle in the current object
              * @return hbond_coefficients_ attribute of the current object of this class
              */
            std::vector<double> GetHbondCoefficients();
/**@*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the list of atom types of the current object
              * Set the types_ attribute of the current bond
              * @param types A list of types of atoms of the current object
              */
            void SetTypes(const std::vector<std::string> types);
            /*! \fn
              * A mutator function in order to set the force constant of the current object
              * Set the force_constant_ attribute of the current angle
              * @param force_constant The force constant of the current object
              */
            void SetForceConstant(double force_constant);
            /*! \fn
              * A mutator function in order to set the length of the current object
              * Set the length_ attribute of the current bond
              * @param length The angle of the current object
              */
            void SetLength(double length);
            /*! \fn
              * A mutator function in order to set the description of the current object
              * Set the dscr_ attribute of the current bond
              * @param dscr The description attribute of the current object
              */
            void SetDscr(const std::string dscr);
            /*! \fn
              * A mutator function in order to set the list of hbond coefficients of the current object
              * Set the hbond_coefficients_ attribute of the current bond
              * @param hbond_coefficients A hbond_coefficients of types of atoms of the current object
              */
            void SetHbondCoefficients(std::vector<double> hbond_coefficients);
/**@*/
            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

         private:
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
