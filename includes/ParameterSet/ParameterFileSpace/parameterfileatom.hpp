#ifndef PARAMETERFILEATOMTYPE_HPP
#define PARAMETERFILEATOMTYPE_HPP

#include <string>
#include <vector>
#include <iostream>

namespace ParameterFileSpace
{
    class ParameterFileAtom
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            ParameterFileAtom();
            /*! \fn
              * Constructor with essential parameters
              * @param type Atom type as a string
              * @param mass A double value of mass of the atom
              * @param polarizability A double value of polarizability of the atom
              * @param dscr A description of the atom that is mentioned in the parameter file
              */
            ParameterFileAtom(const std::string& type, double mass, double polarizability, const std::string& dscr = "");
            /*! \fn
              * Constructor with additional parameters
              * @param type Atom type as a string
              * @param mass Value of mass of the atom
              * @param polarizability Value of polarizability of the atom
              * @param radius Radius of the current atom
              * @param well_depth
              * @param dscr A description of the atom that is mentioned in the parameter file
              * @param mod4_dscr Another description that is mentioned in the parameter file in later section
              */
            ParameterFileAtom(const std::string& type, double mass, double polarizability, double radius,
                              double well_depth, const std::string& dscr = "", const std::string& mod4_dscr = "");           
            /*! \fn
              * Constructor with additional parameters
              * @param type Atom type as a string
              * @param mass Value of mass of the atom
              * @param polarizability Value of polarizability of the atom
              * @param radius Radius of the current atom
              * @param well_depth
              * @param equivalent_list A list of atom types that are equivalent to each other
              * @param dscr A description of the atom that is mentioned in the parameter file
              * @param mod4_dscr Another description that is mentioned in the parameter file in later section
              * @param is_hydrophilic Boolean value to distinguish between hydrophilic and non-hydrophilic atoms
              */
            ParameterFileAtom(const std::string& type, double mass, double polarizability, double radius,
                              double well_depth, const std::vector<std::string>& equivalent_list,
                              const std::string& dscr = "", const std::string& mod4_dscr = "", bool is_hydrophilic = false);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** @addtogroup Molecular_Data_Structure
*  @{
*/
            /*! \fn
              * An accessor function in order to access to type of the current object
              * @return type_ attribute of the current object of this class
              */
            std::string GetType();
            /*! \fn
              * An accessor function in order to access to access to the mass attribute of the current object
              * The attribute is set by the contents of the given file
              * @return mass_ of the current object of this class
              */
            double GetMass();
            /*! \fn
              * An accessor function in order to access to polarizability attribute of the current object
              * The attribute is set by the contents of the given file
              * @return polarizability_ of the current object of this class
              */
            double GetPolarizability();
            /*! \fn
              * An accessor function in order to access to dscr of the current object
              * @return dscr_ attribute of the current object of this class
              */
            std::string GetDscr();
            /*! \fn
              * An accessor function in order to access to radius attribute of the current object
              * The attribute is set by the contents of the given file
              * @return radius_ of the current object of this class
              */
            double GetRadius();
            /*! \fn
              * An accessor function in order to access to well depth attribute of the current object
              * The attribute is set by the contents of the given file
              * @return well_depth_ of the current object of this class
              */
            double GetWellDepth();
            /*! \fn
              * An accessor function in order to access to mod4_dscr of the current object
              * @return mod4_dscr_ attribute of the current object of this class
              */
            std::string GetMod4Dscr();
            /*! \fn
              * An accessor function in order to access to is_hydrophilic attribute of the current object
              * The attribute is set by the contents of the given file
              * @return is_hydrophilic_ of the current object of this class
              */
            bool GetIsHydrophilic();
            /*! \fn
              * An accessor function in order to access to the equivalent list in the current object
              * @return equivalent_list_ attribute of the current object of this class
              */
            std::vector<std::string> GetEquivalentList();
/** @}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** @addtogroup Manipulators
*  @{
*/
            /*! \fn
              * A mutator function in order to set the type of the current object
              * Set the type_ attribute of the current atom
              * @param type The type attribute of the current object
              */
            void SetType( std::string type);
            /*! \fn
              * A mutator function in order to set the mass of the current object
              * Set the mass_ attribute of the current atom
              * @param mass The force mass of the current object
              */
            void SetMass(double mass);
            /*! \fn
              * A mutator function in order to set the polarizability of the current object
              * Set the polarizability_ attribute of the current atom
              * @param polarizability The polarizability of the current object
              */
            void SetPolarizability(double polarizability);
            /*! \fn
              * A mutator function in order to set the description of the current object
              * Set the dscr_ attribute of the current atom
              * @param dscr The description attribute of the current object
              */
            void SetDscr( std::string dscr);
            /*! \fn
              * A mutator function in order to set the radius of the current object
              * Set the radius_ attribute of the current atom
              * @param radius The radius of the current object
              */
            void SetRadius(double radius);
            /*! \fn
              * A mutator function in order to set the well depth of the current object
              * Set the well_depth_ attribute of the current atom
              * @param well_depth The well_depth of the current object
              */
            void SetWellDepth(double well_depth);
            /*! \fn
              * A mutator function in order to set the mod4 dscription of the current object
              * Set the mod4_dscr_ attribute of the current atom
              * @param mod4_dscr The mod4_dscr attribute of the current object
              */
            void SetMod4Dscr( std::string mod4_dscr);
            /*! \fn
              * A mutator function in order to set the is_hydrophilic_ attribute of the current object
              * Set the is_hydrophilic_ attribute of the current atom
              * @param is_hydrophilic The well_depth of the current object
              */
            void SetIsHydrophilic(bool is_hydrophilic);
            /*! \fn
              * A mutator function in order to set equivalent_list of atoms of the current object
              * Set the equivalent_list_ attribute of the current atom
              * @param equivalent_list A list of equivalent_list of atoms of the current object
              */
            void SetEquivalentList( std::vector<std::string> equivalent_list);
/**@}*/
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
            std::string type_;                  /*!< Atom type; Fill by the first column of the first section of the parameter file*/
            double mass_;                       /*!< Atom type mass; Fill by the second column of the first section of the parameter file*/
            double polarizability_;             /*!< Atom type polarizability; Fill by the third column of the first section of the parameter file*/
            std::string dscr_;                  /*!< Atom type description; Fill by the fourth column of the first section of the parameter file*/
            /*!< A sample of the first (Atom Type) section of parameter file: H  1.008         0.161               H bonded to nitrogen atoms*/

            double radius_;                     /*!< Atom type radius; Fill by the corresponding second column of the MOD4 section of the parameter file*/
            double well_depth_;                 /*!< Atom type radius; Fill by the corresponding third column of the MOD4 section of the parameter file*/
            std::string mod4_dscr_;             /*!< Atom type MOD4 description; Fill by the corresponding forth column of the MOD4 section of the parameter file*/
            /*!< A sample of the MOD4 section of the parameter file : H           0.6000  0.0157            !Ferguson base pair geom.*/

            bool is_hydrophilic_;               /*!< Determine hydrophilic atom types; Fill by the second section of the parameter file split by ' '*/
            /*!< A sample of the second (hydrophilic) section of the parameter file: C   H   HO  N   NA  NB  NC  N2  NT  N2  N3  N*  O   OH  OS  P   O2*/

            std::vector<std::string> equivalent_list_;   /*!< Atom types equivalent lists: Fill by the 8th section of the parameter file split by ' '*/
            /*!< A sample of the 8th (equivalent symbols) section of the parameter file: N   NA  N2  N*  NC  NB  NT  NY*/

    };
}

#endif // PARAMETERFILEATOMTYPE_HPP
