#include <sstream>
#include <iomanip>
#include <ios>

#include "includes/common.hpp"
#include "includes/ParameterSet/PrepFile/prepAtom.hpp"

using prep::PrepAtom;
//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepAtom::PrepAtom(std::string& line)
{
	std::string temp;
	std::stringstream ss(line);
	ss >> index_
	>> temp
	>> type_;
	this->setName(temp);
	this->SetTopologicalType(this->ExtractAtomTopologicalType(ss));
	ss >> bond_index_
	>> angle_index_
	>> dihedral_index_
	>> bond_length_
	>> angle_
	>> dihedral_
	>> charge_;
}
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
int PrepAtom::GetIndex() const
{
	return index_;
}

std::string PrepAtom::GetType() const
{
	return type_;
}

gmml::TopologicalType PrepAtom::GetTopologicalType() const
{
	return topological_type_;
}

int PrepAtom::GetBondIndex() const
{
	return bond_index_;
}

int PrepAtom::GetAngleIndex() const
{
	return angle_index_;
}

int PrepAtom::GetDihedralIndex() const
{
	return dihedral_index_;
}

double PrepAtom::GetBondLength() const
{
	return bond_length_;
}

double PrepAtom::GetAngle() const
{
	return angle_;
}

double PrepAtom::GetDihedral() const
{
	return dihedral_;
}

double PrepAtom::GetCharge() const
{
	return charge_;
}
//std::string PrepAtom::GetStringFormatOfTopologicalType(gmml::TopologicalType topological_type)
//{
//    switch(topological_type)
//    {
//        case gmml::kTopTypeE:
//            return "E";
//        case gmml::kTopTypeS:
//            return "S";
//        case gmml::kTopTypeB:
//            return "B";
//        case gmml::kTopType3:
//            return "3";
//        case gmml::kTopType4:
//            return "4";
//        case gmml::kTopTypeM:
//            return "M";
//        default:
//            return "E";
//    }
//}
std::string PrepAtom::GetStringFormatOfTopologicalType() const
{
	switch(this->GetTopologicalType())
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

gmml::TopologicalType PrepAtom::GetTopologicalTypeFromString(std::string topological_type) const
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
void PrepAtom::SetIndex(int index)
{
	index_ = index;
}

void PrepAtom::SetType(const std::string type)
{
	type_ = type;
}

void PrepAtom::SetTopologicalType(gmml::TopologicalType topological_type)
{
	topological_type_ = topological_type;
}

void PrepAtom::SetBondIndex(int bond_index){
	bond_index_ = bond_index;
}

void PrepAtom::SetAngleIndex(int angle_index){
	angle_index_ = angle_index;
}

void PrepAtom::SetDihedralIndex(int dihedral_index){
	dihedral_index_ = dihedral_index;
}

void PrepAtom::SetBondLength(double bond_length){
	bond_length_ = bond_length;
}

void PrepAtom::SetAngle(double angle){
	angle_ = angle;
}

void PrepAtom::SetDihedral(double dihedral){
	dihedral_ = dihedral;
}

void PrepAtom::SetCharge(double charge){
	charge_ = charge;
}
//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
gmml::TopologicalType PrepAtom::ExtractAtomTopologicalType(std::istream &ss)
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
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepAtom::Print(std::ostream &out)
{
	out << std::setw(3) << this->GetIndex()
        				<< std::setw(6) << this->getName()
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

void PrepAtom::Write(std::ostream &stream) const
{
	stream << std::right << std::setw(2) << this->GetIndex() << " " << std::left << std::setw(4) << this->getName() << " " << std::left << std::setw(3) << this->GetType() << " " << std::setw(1) << this->GetStringFormatOfTopologicalType() << " " << std::right << std::setw(2) << this->GetBondIndex() << " " << std::right << std::setw(2) << this->GetAngleIndex() << " " << std::right << std::setw(2) << this->GetDihedralIndex() << " ";
	stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->GetBondLength() << " ";
	stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->GetAngle() << " ";
	stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->GetDihedral();
	stream << "    " << std::right << std::setw(8) << std::fixed << std::setprecision(4) << this->GetCharge() << std::endl;
}

