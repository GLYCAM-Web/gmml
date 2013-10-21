#ifndef PREPFILEATOM_HPP
#define PREPFILEATOM_HPP

#include <string>
#include <iostream>
#include <iostream>

namespace PrepFileSpace
{
    /*! \enum
      * Topological type enumerator
      */
    enum TopologicalType
    {
        kTopTypeE,
        kTopTypeS,
        kTopTypeB,
        kTopType3,
        kTopType4,
        kTopTypeM
    };

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
            PrepFileAtom(int index, const std::string& name, const std::string& type, TopologicalType topological_type, int bond_index,
                         int angle_index, int dihedral_index, double bond_length, double angle, double dihedral, double charge);
            /*! \fn
              * Constructor that extracts information of the atom from the given line from the atom section of each residue in a prep file
              * @param line A line of atom section of each residue in the current prep file
              */
            PrepFileAtom(std::string& line);

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse and extract topological type from a given string stream
              * Parse the given stream and return topological type of the current atom
              * @param ss A string stream that contains the topological type of the atom
              * @return Topological type of the current atom that is in the given string stream
              */
            TopologicalType ExtractAtomTopologicalType(std::istream& ss);

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
            int index_;                                 /*!< Atom index; fill by the first column of the residue section of the file */
            std::string name_;                          /*!< Atom name; fill by the second column of the residue section of the file */
            std::string type_;                          /*!< Atom type; fill by the third column of the residue section of the file */
            TopologicalType topological_type_;          /*!< Topological type (for chain extraction of the residue); fill by th 4th column of the residue section of the file */
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
