#include <sstream>
#include <iomanip>

#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"

using namespace std;
using namespace gmml;
using namespace PrepFileSpace;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepFileAtom::PrepFileAtom() : index_(0), name_(""), type_(""), topological_type_(kTopTypeM), bond_index_(0), angle_index_(0), dihedral_index_(0),
    bond_length_(kNotSet), angle_(kNotSet), dihedral_(kNotSet), charge_(kNotSet) {}

PrepFileAtom::PrepFileAtom(int index, const string& name, const string& type, TopologicalType topological_type, int bond_index, int angle_index, int dihedral_index,
                           double bond_length, double angle, double dihedral, double charge) :
    index_(index), name_(name), type_(type), topological_type_(topological_type), bond_index_(bond_index), angle_index_(angle_index), dihedral_index_(dihedral_index),
    bond_length_(bond_length), angle_(angle), dihedral_(dihedral), charge_(charge) {}

/// Create a prep file atom by a formatted line
PrepFileAtom::PrepFileAtom(std::string& line)
{
    std::stringstream ss(line);
    ss >> index_
       >> name_
       >> type_;

    topological_type_ = ExtractAtomTopologicalType(ss);

    ss >> bond_index_
       >> angle_index_
       >> dihedral_index_
       >> bond_length_
       >> angle_
       >> dihedral_
       >> charge_;
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Extract corresponding topological type from a stream line
/// Return value is an enumeration type
TopologicalType PrepFileAtom::ExtractAtomTopologicalType(istream &ss)
{
    string s;
    ss >> s;
    if (s == "M")
        return kTopTypeM;
    else if (s == "S")
        return kTopTypeS;
    else if (s == "B")
        return kTopTypeB;
    else if (s == "E")
        return kTopTypeE;
    else
        return kTopType3;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepFileAtom::Print(ostream &out)
{
    out << setw(3) << index_
        << setw(6) << name_
        << setw(6) << type_;
    if(topological_type_ == kTopTypeE)
        out << setw(3) << "E";
    else if(topological_type_ == kTopTypeS)
        out << setw(3) << "S";
    else if(topological_type_ == kTopTypeB)
        out << setw(3) << "B";
    else if(topological_type_ == kTopType3)
        out << setw(3) << "3";
    else if(topological_type_ == kTopType4)
        out << setw(3) << "4";
    else if(topological_type_ == kTopTypeM)
        out << setw(3) << "M";
    else
        out << setw(3) << "-";

    out << setw(4) << bond_index_
        << setw(4) << angle_index_
        << setw(4) << dihedral_index_
        << setw(10) << bond_length_
        << setw(10) << angle_
        << setw(10) << dihedral_
        << setw(10) << charge_
        << endl;

}
