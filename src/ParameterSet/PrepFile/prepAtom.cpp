#include <sstream>
#include <iomanip>
#include <ios>

#include "includes/common.hpp" // ToDo find out what's necessary in here and move to dedicated file in CodeUtils
#include "includes/ParameterSet/PrepFile/prepAtom.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp" //get_cartesian_point_from_internal_coords()
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"


using prep::PrepAtom;
//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepAtom::PrepAtom(const std::string& line)
{
	std::stringstream ss(line);
	this->setIndex(codeUtils::extractFromStream(ss, int()));
	this->setName(codeUtils::extractFromStream(ss, std::string()));
	this->SetType(codeUtils::extractFromStream(ss, std::string()));
	this->SetTopologicalType(this->ExtractAtomTopologicalType(ss));
	this->SetBondIndex(codeUtils::extractFromStream(ss, int()));
	this->SetAngleIndex(codeUtils::extractFromStream(ss, int()));
	this->SetDihedralIndex(codeUtils::extractFromStream(ss, int()));
	this->SetBondLength(codeUtils::extractFromStream(ss, double()));
	this->SetAngle(codeUtils::extractFromStream(ss, double()));
	this->SetDihedral(codeUtils::extractFromStream(ss, double()));
	this->setCharge(codeUtils::extractFromStream(ss, double()));
}
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
std::string PrepAtom::GetType() const
{
	return type_;
}

prep::TopologicalType PrepAtom::GetTopologicalType() const
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

std::string PrepAtom::GetStringFormatOfTopologicalType() const
{
	switch(this->GetTopologicalType())
	{
	case kTopTypeE:
		return "E";
	case kTopTypeS:
		return "S";
	case kTopTypeB:
		return "B";
	case kTopType3:
		return "3";
	case kTopType4:
		return "4";
	case kTopTypeM:
		return "M";
	default:
		return "E";
	}
}

prep::TopologicalType PrepAtom::GetTopologicalTypeFromString(std::string topological_type) const
{
	if(topological_type.compare("E") == 0)
		return kTopTypeE;
	if(topological_type.compare("S") == 0)
		return kTopTypeS;
	if(topological_type.compare("B") == 0)
		return kTopTypeB;
	if(topological_type.compare("3") == 0)
		return kTopType3;
	if(topological_type.compare("4") == 0)
		return kTopType4;
	if(topological_type.compare("M") == 0)
		return kTopTypeM;
	else
		return kTopTypeE;
}
//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////


void PrepAtom::SetType(const std::string type)
{
	type_ = type;
}

void PrepAtom::SetTopologicalType(TopologicalType topological_type)
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
//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void PrepAtom::FindDihedralAtoms(std::vector<PrepAtom*>& foundAtoms, int currentDepth, const int& targetDepth)
{
	//std::cout << "Depth is " << currentDepth << " with target " << targetDepth << "\n";
	if(currentDepth == targetDepth)
	{
		return;
	}
	PrepAtom* parent = static_cast<PrepAtom*>(foundAtoms.back()->getParents().front()); // Go up the first parent only. Loops may create another, but they should be ignored.
	foundAtoms.push_back(parent);
	this->FindDihedralAtoms(foundAtoms, ++currentDepth, targetDepth);
	return;
}

void PrepAtom::Determine3dCoordinate()
{
	//std::cout << "Determining 3d Coordinates for " << this->getName() << "\n";
	std::vector<PrepAtom*> foundAtoms;
	foundAtoms.push_back(this);
	this->FindDihedralAtoms(foundAtoms);
	if (foundAtoms.at(3)->getCoordinate() == nullptr)
	{
		std::string message = "This atom has no coordinate: " + foundAtoms.at(3)->getName();
		gmml::log(__LINE__,__FILE__,gmml::ERR, message);
		throw std::runtime_error(message);
	}
	this->setCoordinate(
			GeometryTopology::get_cartesian_point_from_internal_coords(
			foundAtoms.at(3)->getCoordinate(),
			foundAtoms.at(2)->getCoordinate(),
			foundAtoms.at(1)->getCoordinate(),
			this->GetAngle(),
			this->GetDihedral(),
			this->GetBondLength()
			)
	);
	return;
}

prep::TopologicalType PrepAtom::ExtractAtomTopologicalType(std::istream &ss)
{
	std::string s;
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

//void PrepAtom::ExtractIndex(std::istream &ss)
//{
//	int index;
//	ss >> index;
//	this->setIndex(index);
//}
//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepAtom::Print(std::ostream &out) const
{
	out << std::setw(3) << this->getIndex()
        				<< std::setw(6) << this->getName()
						<< std::setw(6) << this->GetType();
	if(this->GetTopologicalType() == kTopTypeE)
		out << std::setw(3) << "E";
	else if(topological_type_ == kTopTypeS)
		out << std::setw(3) << "S";
	else if(topological_type_ == kTopTypeB)
		out << std::setw(3) << "B";
	else if(topological_type_ == kTopType3)
		out << std::setw(3) << "3";
	else if(topological_type_ == kTopType4)
		out << std::setw(3) << "4";
	else if(topological_type_ == kTopTypeM)
		out << std::setw(3) << "M";
	else
		out << std::setw(3) << "-";

	out << std::setw(4) << this->GetBondIndex()
			<< std::setw(4) << this->GetAngleIndex()
			<< std::setw(4) << this->GetDihedralIndex()
			<< std::setw(10) << this->GetBondLength()
			<< std::setw(10) << this->GetAngle()
			<< std::setw(10) << this->GetDihedral()
			<< std::setw(10) << this->getCharge();
	//        << endl;
}

void PrepAtom::Write(std::ostream &stream) const
{
	stream << std::right << std::setw(2) << this->getIndex() << " " << std::left << std::setw(4) << this->getName() << " " << std::left << std::setw(3) << this->GetType() << " " << std::setw(1) << this->GetStringFormatOfTopologicalType() << " " << std::right << std::setw(2) << this->GetBondIndex() << " " << std::right << std::setw(2) << this->GetAngleIndex() << " " << std::right << std::setw(2) << this->GetDihedralIndex() << " ";
	stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->GetBondLength() << " ";
	stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->GetAngle() << " ";
	stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->GetDihedral();
	stream << "    " << std::right << std::setw(8) << std::fixed << std::setprecision(4) << this->getCharge() << std::endl;
}

