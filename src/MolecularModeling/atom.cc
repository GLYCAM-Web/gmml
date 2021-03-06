#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/quantommechanicatom.hpp"
#include "../../includes/MolecularModeling/moleculardynamicatom.hpp"
#include "../../includes/MolecularModeling/dockingatom.hpp"
#include "../../includes/MolecularModeling/oligosaccharidedetectionatom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/GeometryTopology/coordinate.hpp"
#include "cmath"

using MolecularModeling::Atom;
//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Atom::Atom()
{
	this->index_ = this->generateAtomIndex();
	// Call to private helper function.
	this->SetAttributes(NULL, "", GeometryTopology::Coordinate::CoordinateVector(), "", "", "", NULL, "", false, "");
	this->SetBFactor(0);
} // end Default Constructor

Atom::Atom(MolecularModeling::Residue* residue, std::string name, GeometryTopology::Coordinate::CoordinateVector coordinates)
{
	this->index_ = this->generateAtomIndex();
	std::stringstream ss;
	ss << name << "_" << this->GetIndex() << "_" << residue->GetName() << "_?_1_?_?_1";
	// Call to private helper function.
	this->SetAttributes(residue, name, coordinates, "", "", "", NULL, ss.str(), false, "");
	//this->SetBFactor(residue[atom]->GetBFactor())
} // end Constructor

Atom::Atom(MolecularModeling::Residue* residue, std::string name, GeometryTopology::Coordinate coordinate)
{
	this->index_ = this->generateAtomIndex();
	std::stringstream ss;
	ss << name << "_" << this->GetIndex() << "_" << residue->GetName() << "_?_1_?_?_1";
	GeometryTopology::Coordinate::CoordinateVector coordinates;
	coordinates.push_back(new GeometryTopology::Coordinate(coordinate.GetX(), coordinate.GetY(), coordinate.GetZ()));
	this->SetAttributes(residue, name, coordinates, "", "", "", NULL, ss.str(), false, "");
} // end Constructor

Atom::Atom(const Atom* atom)
{
	this->index_ = this->generateAtomIndex();
	this->Copy(atom);
} // end Copy Constructor(*)


/*! \todo  Figure out why the constructor below gets the error below and fix it.
 
src/MolecularModeling/atom.cc: In copy constructor 'MolecularModeling::Atom::Atom(const MolecularModeling::Atom&)':
src/MolecularModeling/atom.cc:51:1: warning: base class 'class MolecularModeling::MolecularDynamicAtom' should be explicitly initialized in the copy constructor [-Wextra]
 Atom::Atom(const Atom& atom)
 ^
*/
// Added by Ayush
Atom::Atom(const Atom& atom)
{
	this->index_ = this->generateAtomIndex();
	this->Copy(&atom);
} // end Copy Constructor(&)

//////////////////////////////////////////////////////////
//                       DESTRUCTOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                         ACCESSORS                    //
//////////////////////////////////////////////////////////
MolecularModeling::Residue* Atom::GetResidue() const
{
	return this->residue_;
} // end GetResidue

std::string Atom::GetName() const
{
	return this->name_;
} // end GetName

GeometryTopology::Coordinate::CoordinateVector Atom::GetCoordinates() const
{
	return this->coordinates_;
} // end GetCoordinates

// Most of the time, you just want the first coordinate
GeometryTopology::Coordinate* Atom::GetCoordinate()
{
    return this->coordinates_.at(0);
} // end GetCoordinate


std::string Atom::GetChemicalType() const
{
 	return this->chemical_type_;
} // end GetChemicalType

std::string Atom::GetDescription() const
{
	return this->description_;
} // end GetDescription

std::string Atom::GetElementSymbol() const
{
	return this->element_symbol_;
} // end GetElementSymbol

MolecularModeling::AtomNode* Atom::GetNode() const
{
	return this->node_;
} // end GetNode

std::string Atom::GetId() const
{
	return this->id_;
} // end GetId

bool Atom::GetIsRing() const
{
	return this->is_ring_;
} // end GetIsRing

unsigned long long Atom::GetIndex() const
{
	return this->index_;
} // end GetIndex

//Added by ayush on 13/11/17 for molecules in assembly to set the atom type as an attribute like O,H, etc.

std::string Atom::GetAtomType() const
{
	return this->atom_type_;
}

//Added by Dave on 03/23/18 for adding B Factor to ontology

float Atom::GetBFactor() const
{
	return b_factor_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Atom::SetResidue(MolecularModeling::Residue* residue)
{
	this->residue_ = residue;
} // end SetResidue

void Atom::SetName(std::string name)
{
	this->name_ = name;
} // end SetName

void Atom::SetCoordinates(GeometryTopology::Coordinate::CoordinateVector coordinates)
{
	// First need to delete any previous Coordinates, so we don't have any memory leaks.
	for(GeometryTopology::Coordinate::CoordinateVector::iterator it = this->coordinates_.begin(); it != this->coordinates_.end(); it++ )
	{
		GeometryTopology::Coordinate* coordinate = (*it);
		if(coordinate != NULL)
		{
			delete coordinate;
			coordinate = NULL;
		}
	}
	this->coordinates_.clear();
	for(GeometryTopology::Coordinate::CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++) {
		this->coordinates_.push_back(*it);
	}
}// end SetCoordinates

void Atom::AddCoordinate(GeometryTopology::Coordinate* coordinate)
{
	this->coordinates_.push_back(coordinate);
} // end AddCoordinate

void Atom::SetChemicalType(std::string chemical_type)
{
	this->chemical_type_ = chemical_type;
} // end SetChemicalType

void Atom::SetDescription(std::string description)
{
	this->description_ = description;
} // end SetDescription

void Atom::SetElementSymbol(std::string element_symbol)
{
	this->element_symbol_ = element_symbol;
} // end SetElementSymbol

void Atom::SetNode(MolecularModeling::AtomNode* node)
{
	this->node_ = node;
} // end SetNode

void Atom::SetId(std::string id)
{
	this->id_ = id;
} // end SetId

void Atom::SetIsRing(bool is_ring)
{
	this->is_ring_ = is_ring;
}
//Added by ayush on 13/11/17 for molecules in assembly

void Atom::SetAtomType(std::string atom_type)
{
	this->atom_type_ = atom_type;
}

//Added by Dave on 03/23/18 for adding B Factor to ontology

void Atom::SetBFactor(float b_factor)
{
	this->b_factor_ = b_factor;
}


unsigned long long Atom::generateAtomIndex() {
	static unsigned long long s_AtomIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
	return s_AtomIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the value in the copy
} // end generateAtomIndex

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Atom::FindConnectedAtoms(AtomVector& visitedAtoms)
{
	visitedAtoms.push_back(this);
	AtomVector neighbors = this->GetNode()->GetNodeNeighbors();
	bool alreadyVisited = false;
	for(AtomVector::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); neighbor++)
	{
		alreadyVisited = false; // reset for each neighbor
		for(AtomVector::iterator visitedAtom = visitedAtoms.begin(); visitedAtom != visitedAtoms.end(); visitedAtom++)
		{
			if((*neighbor)->GetIndex() == (*visitedAtom)->GetIndex())
			alreadyVisited = true;
		}
		if(!alreadyVisited)
		{
            (*neighbor)->FindConnectedAtoms(visitedAtoms); // recursive function call
		}
	}
} // end FindConnectedAtoms

double Atom::GetDistanceToAtom(Atom* otherAtom)
{
    return GetDistanceToCoordinate(otherAtom->GetCoordinate());
} // end GetDistanceToAtom

double Atom::GetDistanceToCoordinate(GeometryTopology::Coordinate* coordinate)
{
    double x = (this->GetCoordinate()->GetX() - coordinate->GetX());
    double y = (this->GetCoordinate()->GetY() - coordinate->GetY());
    double z = (this->GetCoordinate()->GetZ() - coordinate->GetZ());
	return sqrt((x * x) + (y * y) + (z * z));
} // end GetDistanceToCoordinate




//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Atom::Print(std::ostream& out)
{
	out << "Atom name: "    << this->GetName() << std::endl;
	out << "Element: "      << this->GetElementSymbol() << std::endl;
	out << "Atom Type: "    << this->GetAtomType() << std::endl;
	out << "Coordinates: "  << std::endl;
	GeometryTopology::Coordinate::CoordinateVector coordinates = this->GetCoordinates();
	for(GeometryTopology::Coordinate::CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
	{
		GeometryTopology::Coordinate* coordinate = (*it);
		if(coordinate != NULL)
		{
			out << "\t";
			coordinate->Print(out);
			out << std::endl;
		}
	}
	out << "**************** Structure *****************" << std::endl;
	if(this->GetNode() != NULL)
	{
		this->GetNode()->Print(out);
	}
	out << std::endl;
} // end Print

//////////////////////////////////////////////////////////
//                   HELPER FUNCTIONS                   //
//////////////////////////////////////////////////////////
void Atom::Copy(const Atom* atom)
{
	// Copy the easy stuff.
	this->SetName(atom->GetName());
	this->SetChemicalType(atom->GetChemicalType());
	this->SetDescription(atom->GetDescription());
	this->SetElementSymbol(atom->GetElementSymbol());
	this->SetId(atom->GetId());
	this->SetIsRing(atom->GetIsRing());
	this->SetAtomType(atom->GetAtomType());
	// Deep Copy objects
	// this->residue_ = new MolecularModeling::Residue(atom->GetResidue());
	this->SetResidue(atom->GetResidue());
	if(this->node_ != NULL)
	{
		delete this->node_;
	}
	this->node_ = new MolecularModeling::AtomNode(atom->GetNode());
	GeometryTopology::Coordinate::CoordinateVector atomCoordinates = atom->GetCoordinates();
	for(GeometryTopology::Coordinate::CoordinateVector::iterator it = atomCoordinates.begin(); it != atomCoordinates.end(); it++ )
	{
		GeometryTopology::Coordinate* tempCoordinate = (*it);
		this->coordinates_.push_back(new GeometryTopology::Coordinate(tempCoordinate));
	}
} // end Copy

void Atom::SetAttributes(	MolecularModeling::Residue* residue, std::string name, GeometryTopology::Coordinate::CoordinateVector coordinates,
							std::string chemical_type, std::string description, std::string element_symbol,
							MolecularModeling::AtomNode* node, std::string id, bool is_ring, std::string atom_type)
{
	// Having this function call the Setter functions for everything allows for
	//	simple error handling and debugging because the setting of the variables
	//	only happen in one place. DT
	this->SetResidue(residue);
	this->SetName(name);
	this->SetCoordinates(coordinates);
	this->SetChemicalType(chemical_type);
	this->SetDescription(description);
	this->SetElementSymbol(element_symbol);
	this->SetNode(node);
	this->SetId(id);
	this->SetIsRing(is_ring);
	// This function doesn't set index because of the attributes uniqueness, it should
	//	only be called from the Constructors. DT
	this->SetAtomType(atom_type);
} // end SetAttributes
// @TODO DT - The below operators are left for future development. See atom.hpp for
//			details on why they could be useful.
// void Atom::operator=(const Atom&)
// {
//
// } // end Atom::operator=(&)

// void Atom::operator=(const Atom*)
// {
//
// } // end operator=(*)

// std::ostream& operator<<(std::ostream& out, const Atom& atom)
// {
// 	out << "Atom name: "    << atom.GetName() << std::endl;
// 	out << "Element: "      << atom.GetElementSymbol() << std::endl;
// 	out << "Atom Type: "    << atom.MolecularDynamicAtom::GetAtomType()	<< std::endl;
// 	out << "Coordinates: "  << std::endl;
// 	for( CoordinateVector::iterator it = atom.GetCoordinates().begin(); it != atom.GetCoordinates().end(); it++ )
// 	{
// 		GeometryTopology::Coordinate * coordinate = (*it);
// 		out << "\t";
// 		coordinate->Print(out);
// 		out << std::endl;
// 	}
// 	if(atom.GetNode() != NULL)
// 	{
// 		out << "**************** Structure *****************" << std::endl;
// 		atom.GetNode()->Print(out);
// 	}
// 	return out;
// } // end opeartor<<
