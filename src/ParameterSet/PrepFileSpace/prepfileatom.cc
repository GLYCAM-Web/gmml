#include <sstream>
#include <iomanip>

#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"

using PrepFileSpace::PrepFileAtom;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepFileAtom::PrepFileAtom() : index_(0), name_(""), type_(""), topological_type_(gmml::kTopTypeM), bond_index_(0), angle_index_(0), dihedral_index_(0),
    bond_length_(gmml::dNotSet), angle_(gmml::dNotSet), dihedral_(gmml::dNotSet), charge_(gmml::dNotSet) {}

PrepFileAtom::PrepFileAtom(int index, const std::string& name, const std::string& type, gmml::TopologicalType topological_type, int bond_index, int angle_index, int dihedral_index,
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

PrepFileAtom::~PrepFileAtom()
{

}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Extract corresponding topological type from a stream line
/// Return value is an enumeration type
gmml::TopologicalType PrepFileAtom::ExtractAtomTopologicalType(std::istream &ss)
{
    std::string s;
    ss >> s;
    if (s == "M")
        return gmml::kTopTypeM;
    else if (s == "S")
        return gmml::kTopTypeS;
    else if (s == "B")
        return gmml::kTopTypeB;
    else if (s == "E")
        return gmml::kTopTypeE;
    else
        return gmml::kTopType3;
}
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////

int PrepFileAtom::GetIndex(){
    return index_;
}

std::string PrepFileAtom::GetName(){
    return name_;
}

std::string PrepFileAtom::GetType(){
    return type_;
}

gmml::TopologicalType PrepFileAtom::GetTopologicalType(){
    return topological_type_;
}

int PrepFileAtom::GetBondIndex(){
    return bond_index_;
}

int PrepFileAtom::GetAngleIndex(){
    return angle_index_;
}

int PrepFileAtom::GetDihedralIndex(){
    return dihedral_index_;
}

double PrepFileAtom::GetBondLength(){
    return bond_length_;
}

double PrepFileAtom::GetAngle(){
    return angle_;
}

double PrepFileAtom::GetDihedral(){
    return dihedral_;
}

double PrepFileAtom::GetCharge(){
    return charge_;
}
std::string PrepFileAtom::GetStringFormatOfTopologicalType(gmml::TopologicalType topological_type)
{
    switch(topological_type)
    {
        case gmml::kTopTypeE:
            return "E";
        case gmml::kTopTypeS:
            return "S";
        case gmml::kTopTypeB:
            return "B";
        case gmml::kTopType3:
            return "3";
        case gmml::kTopType4:
            return "4";
        case gmml::kTopTypeM:
            return "M";
        default:
            return "E";
    }
}
std::string PrepFileAtom::GetStringFormatOfTopologicalType()
{
    switch(topological_type_)
    {
        case gmml::kTopTypeE:
            return "E";
        case gmml::kTopTypeS:
            return "S";
        case gmml::kTopTypeB:
            return "B";
        case gmml::kTopType3:
            return "3";
        case gmml::kTopType4:
            return "4";
        case gmml::kTopTypeM:
            return "M";
        default:
            return "E";
    }
}
gmml::TopologicalType PrepFileAtom::GetTopologicalTypeFromString(std::string topological_type)
{
    if(topological_type.compare("E") == 0)
        return gmml::kTopTypeE;
    if(topological_type.compare("S") == 0)
        return gmml::kTopTypeS;
    if(topological_type.compare("B") == 0)
        return gmml::kTopTypeB;
    if(topological_type.compare("3") == 0)
        return gmml::kTopType3;
    if(topological_type.compare("4") == 0)
        return gmml::kTopType4;
    if(topological_type.compare("M") == 0)
        return gmml::kTopTypeM;
    else
        return gmml::kTopTypeE;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////

void PrepFileAtom::SetIndex(int index){
    index_ = index;
}

void PrepFileAtom::SetName(const std::string name){
    name_ = name;
}

void PrepFileAtom::SetType(const std::string type){
    type_ = type;
}

void PrepFileAtom::SetTopologicalType(gmml::TopologicalType topological_type){
    topological_type_ = topological_type;
}

void PrepFileAtom::SetBondIndex(int bond_index){
    bond_index_ = bond_index;
}

void PrepFileAtom::SetAngleIndex(int angle_index){
    angle_index_ = angle_index;
}

void PrepFileAtom::SetDihedralIndex(int dihedral_index){
    dihedral_index_ = dihedral_index;
}

void PrepFileAtom::SetBondLength(double bond_length){
    bond_length_ = bond_length;
}

void PrepFileAtom::SetAngle(double angle){
    angle_ = angle;
}

void PrepFileAtom::SetDihedral(double dihedral){
    dihedral_ = dihedral;
}

void PrepFileAtom::SetCharge(double charge){
    charge_ = charge;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepFileAtom::Print(std::ostream &out)
{
    out << std::setw(3) << index_
        << std::setw(6) << name_
        << std::setw(6) << type_;
    if(topological_type_ == gmml::kTopTypeE)
        out << std::setw(3) << "E";
    else if(topological_type_ == gmml::kTopTypeS)
        out << std::setw(3) << "S";
    else if(topological_type_ == gmml::kTopTypeB)
        out << std::setw(3) << "B";
    else if(topological_type_ == gmml::kTopType3)
        out << std::setw(3) << "3";
    else if(topological_type_ == gmml::kTopType4)
        out << std::setw(3) << "4";
    else if(topological_type_ == gmml::kTopTypeM)
        out << std::setw(3) << "M";
    else
        out << std::setw(3) << "-";

    out << std::setw(4) << bond_index_
        << std::setw(4) << angle_index_
        << std::setw(4) << dihedral_index_
        << std::setw(10) << bond_length_
        << std::setw(10) << angle_
        << std::setw(10) << dihedral_
        << std::setw(10) << charge_;
//        << endl;

}
