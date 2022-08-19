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

double PrepAtom::GetCharge() const
{
	return charge_;
}
//std::string PrepAtom::GetStringFormatOfTopologicalType(TopologicalType topological_type)
//{
//    switch(topological_type)
//    {
//        case kTopTypeE:
//            return "E";
//        case kTopTypeS:
//            return "S";
//        case kTopTypeB:
//            return "B";
//        case kTopType3:
//            return "3";
//        case kTopType4:
//            return "4";
//        case kTopTypeM:
//            return "M";
//        default:
//            return "E";
//    }
//}
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
void PrepAtom::SetIndex(int index)
{
	index_ = index;
}

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

void PrepAtom::SetCharge(double charge){
	charge_ = charge;
}
//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void PrepAtom::FindDihedralAtoms(std::vector<PrepAtom*>& foundAtoms, int currentDepth, const int& targetDepth)
{


//	 ToDo PrepAtom::FindDihedralAtoms is the function that exposed the weakness in the design. I can't use Node class functions here, I get back
//	 the wrong type ie:
//	 error: invalid conversion from ‘__gnu_cxx::__alloc_traits<std::allocator<glygraph::Node<cds::cdsAtom>*>, glygraph::Node<cds::cdsAtom>*>::value_type’ {aka ‘glygraph::Node<cds::cdsAtom>*’} to ‘prep::PrepAtom*’
//	 There are two issues. getParents should return:
//	std::vector<T*> Node<T>::getParents() and not std::vector<Node<T> *> Node<T>::getParents()
//	 The other issue is that prepAtom inherits from glygraph::Node<cds::cdsAtom>, but it should instead be a node ie.
//	 PrepAtom : public glygraph::Node<prepAtom>, which is the CRTP we are using everywhere. It can get atom attributes
//	 the same as we give Residue attributes ie. PrepAtom : public glygraph::Node<prepAtom>, public Abstract::Atom
//	 I can convert cds::Atom into Abstrac::Atom and inherit it in pdbAtom and prepAtom.
//	 Another option is to have cds::Atom be a template like cds::Residue and higher classes are, but I can't conceptualize it
//	 properly. cdsResidue contains atoms whose type is templated. So it cna have prepAtoms or pdbAtoms or whatever when
//	 it becomes a prepResidue like here: class PrepResidue : public cds::cdsResidue<PrepAtom>
//     but I don't get what I should have cdsAtom be. Coordinate is concrete so there's nothing to template. I what I need is
//	 for it to be Abstract::Atom...




	//PrepAtom* firstParent = foundAtoms.back()->getParents().front(); // Go up the first one only
	//std::cout << firstParent->getName() << std::endl;
	//std::cout << firstParent->GetType() << std::endl;
//	foundAtoms.push_back(firstParent);
//	if(currentDepth == targetDepth)
//	{
//		return;
//	}
//	this->FindDihedralAtoms(foundAtoms, ++currentDepth, targetDepth);
//	return;
}


void PrepAtom::Determine3dCoordinate()
{
	std::cout << "Determining 3d Coordinates!\n";
	std::vector<PrepAtom*> foundAtoms;
	foundAtoms.push_back(this);
	this->FindDihedralAtoms(foundAtoms);
	std::cout << "Found these atoms:\n";
	for(auto &dihedralAtom : foundAtoms)
	{
		std::cout << dihedralAtom->getName() << std::endl;
	}
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
//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepAtom::Print(std::ostream &out) const
{
	out << std::setw(3) << this->GetIndex()
        				<< std::setw(6) << this->getName()
						<< std::setw(6) << type_;
	if(topological_type_ == kTopTypeE)
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

