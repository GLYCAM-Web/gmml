#ifndef PREPFILEATOM_HPP
#define PREPFILEATOM_HPP

#include <string>
#include <iostream>
#include <iostream>
#include "../../common.hpp"

namespace PrepFileSpace
{
    class PrepFileAtom; // Forward declare for the vector typedef:
    typedef std::vector<PrepFileAtom*> PrepFileAtomVector;
    class PrepFileAtom
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PrepFileAtom();
            /*! \fn
              * Constructor with essential parameters
              * @param index Atom index in its belonging residue
              * @param name Atom name
              * @param type Atom type
              * @param topological_type Topological type of the atom in its residue that can be used to extract residue structure
              * @param bond_index Atom index in the bond
              * @param angle_index Atom index in the angle
              * @param dihedral_index Atom index in the dihedral
              * @param bond_length The actual value of the bond that the atom is involved in
              * @param angle The actual value of the angle that the atom is involved in
              * @param dihedral The actual value of the dihedral that the atom is involved in
              * @param charge The actual value of charge of the atom
              */
            PrepFileAtom(int index, const std::string& name, const std::string& type, gmml::TopologicalType topological_type, int bond_index,
                         int angle_index, int dihedral_index, double bond_length, double angle, double dihedral, double charge);
            /*! \fn
              * Constructor that extracts information of the atom from the given line from the atom section of each residue in a prep file
              * @param line A line of atom section of each residue in the current prep file
              */
            PrepFileAtom(std::string& line);

            ~PrepFileAtom();
            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * A function to parse and extract topological type from a given string stream
              * Parse the given stream and return topological type of the current atom
              * @param ss A string stream that contains the topological type of the atom
              * @return Topological type of the current atom that is in the given string stream
              */
            gmml::TopologicalType ExtractAtomTopologicalType(std::istream& ss);
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

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to index of the current object
              * @return index_ attribute of the current object of this class
              */
            int GetIndex();
            /*! \fn
              * An accessor function in order to access to access to the name attribute of the current object
              * The attribute is set by the contents of the given file
              * @return name_ of the current object of this class
              */
            std::string GetName();
            /*! \fn
              * An accessor function in order to access to type attribute of the current object
              * The attribute is set by the contents of the given file
              * @return type_ of the current object of this class
              */
            std::string GetType();
            /*! \fn
              * An accessor function in order to access to topological type of the current object
              * @return topological_type_ attribute of the current object of this class
              */
            gmml::TopologicalType GetTopologicalType();
            /*! \fn
              * An accessor function in order to access to bond index of the current object
              * The attribute is set by the contents of the given file
              * @return bond_index_ of the current object of this class
              */
            int GetBondIndex();
            /*! \fn
              * An accessor function in order to access to angle index of the current object
              * The attribute is set by the contents of the given file
              * @return angle_index_ of the current object of this class
              */
            int GetAngleIndex();
            /*! \fn
              * An accessor function in order to access to dihedral index of the current object
              * The attribute is set by the contents of the given file
              * @return dihedral_index_ of the current object of this class
              */
            int GetDihedralIndex();
            /*! \fn
              * An accessor function in order to access to bond length of the current object
              * The attribute is set by the contents of the given file
              * @return bond_length_ of the current object of this class
              */
            double GetBondLength();
            /*! \fn
              * An accessor function in order to access to angle attribute of the current object
              * The attribute is set by the contents of the given file
              * @return angle_ of the current object of this class
              */
            double GetAngle();
            /*! \fn
              * An accessor function in order to access to dihedral attribute of the current object
              * The attribute is set by the contents of the given file
              * @return dihedral_ of the current object of this class
              */
            double GetDihedral();
            /*! \fn
              * An accessor function in order to access to charge attribute of the current object
              * The attribute is set by the contents of the given file
              * @return charge_ of the current object of this class
              */
            double GetCharge();
            /*! \fn
              * Convert a value of TopologicalType enumerator to the string version of it
              * @param topological_type A value of TopologicalType has to be converted to string
              * @return String format of the given value of TopologicalType enumerator
              */
            std::string GetStringFormatOfTopologicalType(gmml::TopologicalType topological_type);
            /*! \fn
              * Convert the value of TopologicalType attribute of the current object (topological_type_) to the string version of it
              * @return String format of the value of topological_type_ attribute of the current object
              */
            std::string GetStringFormatOfTopologicalType();
            /*! \fn
              * Convert string version of TopologicalType to the corresponding enum value
              * @param topological_type String indicates TopologicalType
              * @return A value selected from TopologicalType enumerator correspondence to the given string
              */
            gmml::TopologicalType GetTopologicalTypeFromString(std::string topological_type);
/**@}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the type of the current object
              * Set the index_ attribute of the current object
              * @param index The index attribute of the current object
              */
            void SetIndex(int index);
            /*! \fn
              * A mutator function in order to set the name of the current object
              * Set the name_ attribute of the current object
              * @param name The force name of the current object
              */
            void SetName(std::string name);
            /*! \fn
              * A mutator function in order to set the type of the current object
              * Set the type_ attribute of the current object
              * @param type The type of the current object
              */
            void SetType(std::string type);
            /*! \fn
              * A mutator function in order to set the topological type of the current object
              * Set the topological_type_ attribute of the current object
              * @param topological_type The topological_type attribute of the current object
              */
            void SetTopologicalType(gmml::TopologicalType topological_type);
            /*! \fn
              * A mutator function in order to set the bond index of the current object
              * Set the bond_index_ attribute of the current object
              * @param bond_index The bond_index of the current object
              */
            void SetBondIndex(int bond_index);
            /*! \fn
              * A mutator function in order to set the angle index of the current object
              * Set the angle_index_ attribute of the current object
              * @param angle_index The angle_index of the current object
              */
            void SetAngleIndex(int angle_index);
            /*! \fn
              * A mutator function in order to set the dihedral index of the current object
              * Set the dihedral_index_ attribute of the current atom
              * @param dihedral_index The dihedral_index attribute of the current object
              */
            void SetDihedralIndex(int dihedral_index);
            /*! \fn
              * A mutator function in order to set the bond length of the current object
              * Set the bond_length_ attribute of the current object
              * @param bond_length The bond_length of the current object
              */
            void SetBondLength(double bond_length);
            /*! \fn
              * A mutator function in order to set the angle of the current object
              * Set the angle_ attribute of the current object
              * @param angle The angle of the current object
              */
            void SetAngle(double angle);
            /*! \fn
              * A mutator function in order to set the dihedral of the current object
              * Set the dihedral_ attribute of the current object
              * @param dihedral The dihedral of the current object
              */
            void SetDihedral(double dihedral);
            /*! \fn
              * A mutator function in order to set the charge of the current object
              * Set the charge_ attribute of the current object
              * @param charge The charge of the current object
              */
            void SetCharge(double charge);
/**@}*/

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            int index_;                                 /*!< Atom index; fill by the first column of the residue section of the file */
            std::string name_;                          /*!< Atom name; fill by the second column of the residue section of the file */
            std::string type_;                          /*!< Atom type; fill by the third column of the residue section of the file */
            gmml::TopologicalType topological_type_;          /*!< Topological type (for chain extraction of the residue); fill by th 4th column of the residue section of the file */
            int bond_index_;                            /*!< Bond index; fill by the 5th column of the residue section of the file */
            int angle_index_;                           /*!< Angle index; fill by the 6th column of the residue section of the file */
            int dihedral_index_;                        /*!< Dihedral index; fill by the 7th column of the residue section of the file */
            double bond_length_;                        /*!< Bond length; fill by the 8th column of the residue section of the file */
            double angle_;                              /*!< Angle; fill by the 9th column of the residue section of the file */
            double dihedral_;                           /*!< Dihedral; fill by the 10th column of the residue section of the file */
            double charge_;                             /*!< Charge; fill by the 11th column of the residue section of the file */
            /*!< Sample line of the residue section of a prep file: 4 H1   H1  M  3  2  1  1.000    90.0     180.0     0.0 */
    };
}

#endif // PREPFILEATOM_HPP
