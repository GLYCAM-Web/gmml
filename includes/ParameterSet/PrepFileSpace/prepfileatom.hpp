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
            PrepFileAtom(int index, const std::string& name, const std::string& type, TopologicalType topological_type, int bond_index, int angle_index, int dihedral_index,
                         double bond_length, double angle, double dihedral, double charge);
            PrepFileAtom(std::string& line);

            ///////////////////////////// FUNCTIONS ///////////////////////////////////
            TopologicalType ExtractAtomTopologicalType(std::istream& ss);

            ////////////////////////// DISPLAY FUNCTION ///////////////////////////////
            void Print(std::ostream& out);

            ///////////////////////////// ATTRIBUTES //////////////////////////////////
            int index_;
            std::string name_;
            std::string type_;
            TopologicalType topological_type_;
            int bond_index_;
            int angle_index_;
            int dihedral_index_;
            double bond_length_;
            double angle_;
            double dihedral_;
            double charge_;
    };
}

#endif // PREPFILEATOM_HPP
