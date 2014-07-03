#include "../../includes/MolecularModeling/atom.hpp"

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
Geometry::Coordinate Atom::GetCoordinate()
{
    return coordinate_;
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
void Atom::SetCoordinate(Geometry::Coordinate coordinate)
{
    coordinate_ = coordinate;
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

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Atom::Print(ostream &out)
{
}
