#ifndef PREPFILEATOM_HPP
#define PREPFILEATOM_HPP

#include <string>

namespace PrepFileSpace
{
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
            /////////////////////////////// CONSTRUCTOR ////////////////////////////////
            PrepFileAtom();
            PrepFileAtom(int index, const std::string& name, const std::string& type, TopologicalType topological_type, int bond_index,
                         int angle_index, int dihedral_index, double bond_length, double angle, double dihedral, double charge);
            PrepFileAtom(std::string& line);

            ///////////////////////////// FUNCTIONS ///////////////////////////////////
            TopologicalType ExtractAtomTopologicalType(std::istream& ss);

            ////////////////////////// DISPLAY FUNCTION ///////////////////////////////
            void Print(std::ostream& out);

            ///////////////////////////// ATTRIBUTES //////////////////////////////////
            int index_;                                 // Atom index; fill by the first column of the residue section of the file
            std::string name_;                          // Atom name; fill by the second column of the residue section of the file
            std::string type_;                          // Atom type; fill by the third column of the residue section of the file
            TopologicalType topological_type_;          // Topological type (for chain extraction of the residue); fill by th 4th column of the residue section of the file
            int bond_index_;                            // Bond index; fill by the 5th column of the residue section of the file
            int angle_index_;                           // Angle index; fill by the 6th column of the residue section of the file
            int dihedral_index_;                        // Dihedral index; fill by the 7th column of the residue section of the file
            double bond_length_;                        // Bond length; fill by the 8th column of the residue section of the file
            double angle_;                              // Angle; fill by the 9th column of the residue section of the file
            double dihedral_;                           // Dihedral; fill by the 10th column of the residue section of the file
            double charge_;                             // Charge; fill by the 11th column of the residue section of the file
            // Sample line of the residue section of a prep file: 4 H1   H1  M  3  2  1  1.000    90.0     180.0     0.0
    };
}

#endif // PREPFILEATOM_HPP
