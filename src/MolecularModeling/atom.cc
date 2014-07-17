#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/quantommechanicatom.hpp"
#include "../../includes/MolecularModeling/moleculardynamicatom.hpp"
#include "../../includes/MolecularModeling/dockingatom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Atom::Atom() {}

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
void Atom::AddCoordinate(Geometry::Coordinate *coordinate)
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

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Atom::Print(ostream &out)
{
}
