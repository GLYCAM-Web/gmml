#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/quantommechanicatom.hpp"
#include "../../includes/MolecularModeling/moleculardynamicatom.hpp"
#include "../../includes/MolecularModeling/dockingatom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "cmath"

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Atom::Atom() : name_(""), chemical_type_(""), element_symbol_(""), description_(""), id_("")
{
    coordinates_ = CoordinateVector();
    residue_ = NULL;
    node_ = NULL;
    index_ = this->generateAtomIndex();
}

Atom::Atom(Residue *residue, string name, CoordinateVector coordinates) :
    chemical_type_(""), element_symbol_(""), description_("")
{
    residue_ = residue;
    name_ = name;
    coordinates_ = CoordinateVector();
    for(CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
        coordinates_.push_back(*it);
    node_ = NULL;
    index_ = this->generateAtomIndex();
}

Atom::Atom(Atom *atom)
{
    residue_ = new Residue(atom->GetResidue());
    name_ = atom->GetName();
    coordinates_ = CoordinateVector();
    CoordinateVector coordinates = atom->GetCoordinates();
    for(CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
        coordinates_.push_back(new GeometryTopology::Coordinate(*it));

    AtomNode node = atom->GetNode();
    node_ = new AtomNode(node);
    index_ = atom->GetIndex();
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Residue* Atom::GetResidue()
{
    return residue_;
}
string Atom::GetName()
{
    return name_;
}
Atom::CoordinateVector Atom::GetCoordinates()
{
    return coordinates_;
}
string Atom::GetChemicalType()
{
    return chemical_type_;
}
string Atom::GetDescription()
{
    return description_;
}
string Atom::GetElementSymbol()
{
    return element_symbol_;
}
AtomNode* Atom::GetNode()
{
    return node_;
}
string Atom::GetId()
{
    return id_;
}
bool Atom::GetIsRing()
{
    return is_ring_;
}
unsigned long long Atom::GetIndex()
{
    return index_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Atom::SetResidue(Residue *residue)
{
    residue_ = residue;
}
void Atom::SetName(string name)
{
    name_ = name;
}
void Atom::SetCoordinates(CoordinateVector coordinates)
{
    coordinates_.clear();
    for(CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
    {
        coordinates_.push_back(*it);
    }
}
void Atom::AddCoordinate(GeometryTopology::Coordinate *coordinate)
{
    coordinates_.push_back(coordinate);
}
void Atom::SetChemicalType(string chemical_type)
{
    chemical_type_ = chemical_type;
}
void Atom::SetDescription(string description)
{
    description_ = description;
}
void Atom::SetElementSymbol(string element_symbol)
{
    element_symbol_ = element_symbol;
}
void Atom::SetNode(AtomNode *node)
{
    node_ = node;
}
void Atom::SetId(string id)
{
    id_ = id;
}
void Atom::SetIsRing(bool is_ring)
{
    is_ring_ = is_ring;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Atom::FindConnectedAtoms(AtomVector &visitedAtoms)
{
    visitedAtoms.push_back(this);
    AtomVector neighbors = this->GetNode()->GetNodeNeighbors();
    bool alreadyVisited = false;

    for(AtomVector::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); neighbor++){
        alreadyVisited = false; // reset for each neighbor
        for(AtomVector::iterator visitedAtom = visitedAtoms.begin(); visitedAtom != visitedAtoms.end(); visitedAtom++){
            if ( (*neighbor)->GetIndex() == (*visitedAtom)->GetIndex() )
                alreadyVisited = true;
        }
        if (!alreadyVisited) {
            //std::cout << "Found unvisited neighbor, Going to " << (*neighbor)->GetId() << " from " << this->GetId() << std::endl;
            (*neighbor)->FindConnectedAtoms(visitedAtoms); // recursive function call
        }
    }
}

double Atom::GetDistanceToAtom(Atom *otherAtom)
{
    double x = ( this->GetCoordinates().at(0)->GetX() - otherAtom->GetCoordinates().at(0)->GetX() );
    double y = ( this->GetCoordinates().at(0)->GetY() - otherAtom->GetCoordinates().at(0)->GetY() );
    double z = ( this->GetCoordinates().at(0)->GetZ() - otherAtom->GetCoordinates().at(0)->GetZ() );
    return sqrt( (x*x) + (y*y) + (z*z) );
}

unsigned long long Atom::generateAtomIndex()
{
    static unsigned long long s_AtomIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
    return s_AtomIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the value in the copy
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Atom::Print(ostream &out)
{
    out << "Atom name: " << name_ << endl;
    out << "Element: " << element_symbol_ << endl;
    out << "Atom Type: " << this->MolecularDynamicAtom::GetAtomType() << endl;
    out << "Coordinates: " << endl;
    if(coordinates_.size() != 0)
    {
        for(CoordinateVector::iterator it = coordinates_.begin(); it != coordinates_.end(); it++)
        {
            GeometryTopology::Coordinate* coordinate = *it;
            out << "\t";
            coordinate->Print(out);
            out << endl;
        }
    }
    out << "**************** Structure *****************" << endl;
    if(node_ != NULL)
        node_->Print(out);
    out << endl;

}
