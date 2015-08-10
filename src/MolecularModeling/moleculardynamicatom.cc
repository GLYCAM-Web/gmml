#include "../../includes/MolecularModeling/moleculardynamicatom.hpp"
#include "../../includes/common.hpp"

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
MolecularDynamicAtom::MolecularDynamicAtom() : atom_type_(""), charge_(gmml::dNotSet), mass_(gmml::dNotSet), radius_(gmml::dNotSet) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string MolecularDynamicAtom::GetAtomType()
{
    return atom_type_;
}
double MolecularDynamicAtom::GetCharge()
{
    return charge_;
}
double MolecularDynamicAtom::GetMass()
{
    return mass_;
}
double MolecularDynamicAtom::GetRadius()
{
    return radius_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void MolecularDynamicAtom::SetAtomType(string atom_type)
{
    atom_type_ = atom_type;
}
void MolecularDynamicAtom::SetCharge(double charge)
{
    charge_ = charge;
}
void MolecularDynamicAtom::SetMass(double mass)
{
    mass_ = mass;
}
void MolecularDynamicAtom::SetRadius(double radius)
{
    radius_ = radius;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void MolecularDynamicAtom::Print(ostream &out)
{    
}

void MolecularDynamicAtom::PrintHet(ostream &out)
{
}
