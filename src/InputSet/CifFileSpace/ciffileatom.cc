
#include "../../../includes/InputSet/CifFileSpace/ciffileatom.hpp"
#include "../../../includes/utils.hpp"

using CifFileSpace::CifFileAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CifFileAtom::CifFileAtom() {}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string CifFileAtom::GetAtomName()
{
    return atom_name_;
}
std::string CifFileAtom::GetAlternateAtomName()
{
    return alternate_atom_name_;
}
std::string CifFileAtom::GetElementSymbol()
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
void CifFileAtom::SetAtomName(std::string atom_name)
{
    atom_name_ = atom_name;
}
void CifFileAtom::SetAlternateAtomName(std::string alternate_atom_name)
{
    alternate_atom_name_ = alternate_atom_name;
}
void CifFileAtom::SetCharge(double charge)
{
    charge_ = charge;
}
void CifFileAtom::SetElementSymbol(std::string element_symbol)
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
void CifFileAtom::Print(std::ostream &out)
{
    out << "";
}
