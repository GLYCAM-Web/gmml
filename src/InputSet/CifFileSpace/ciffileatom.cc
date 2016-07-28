
#include "../../../includes/InputSet/CifFileSpace/ciffileatom.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace CifFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CifFileAtom::CifFileAtom() {}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string CifFileAtom::GetAtomName()
{
    return atom_name_;
}
string CifFileAtom::GetAlternateAtomName()
{
    return alternate_atom_name_;
}
string CifFileAtom::GetElementSymbol()
{
    return element_symbol_;
}

double CifFileAtom::GetCharge()
{
    return charge_;
}

GeometryTopology::Coordinate* CifFileAtom::GetCoordinate()
{
    return coordinate_;
}
GeometryTopology::Coordinate* CifFileAtom::GetCoordinateIdeal()
{
    return coordinate_ideal_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void CifFileAtom::SetAtomName(string atom_name)
{
    atom_name_ = atom_name;
}
void CifFileAtom::SetAlternateAtomName(string alternate_atom_name)
{
    alternate_atom_name_ = alternate_atom_name;
}
void CifFileAtom::SetCharge(double charge)
{
    charge_ = charge;
}
void CifFileAtom::SetElementSymbol(string element_symbol)
{
    element_symbol_ = element_symbol;
}
void CifFileAtom::SetCoordinate(GeometryTopology::Coordinate *coordinate)
{
    coordinate_ = coordinate;
}
void CifFileAtom::SetCoordinateIdeal(GeometryTopology::Coordinate *coordinate_ideal)
{
    coordinate_ideal_ = coordinate_ideal;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void CifFileAtom::Print(ostream &out)
{
}


