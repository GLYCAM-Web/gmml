#ifndef PARAMETERFILEATOMTYPE_HPP
#define PARAMETERFILEATOMTYPE_HPP

#include <string>
#include <vector>

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
