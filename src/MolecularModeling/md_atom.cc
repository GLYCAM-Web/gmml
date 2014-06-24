#include "../../../includes/MolecularModeling/md_atom.hpp"

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
MD_Atom::MD_Atom() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string MD_Atom::GetAtomType()
{
    return atom_type_;
}
double MD_Atom::GetCharge()
{
    return charge_;
}
double MD_Atom::GetMass()
{
    return mass_;
}
double MD_Atom::GetRadius()
{
    return radius_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void MD_Atom::SetAtomType(string atom_type)
{
    atom_type_ = atom_type;
}
void MD_Atom::SetCharge(double charge)
{
    charge_ = charge;
}
void MD_Atom::SetMass(double mass)
{
    mass_ = mass;
}
void MD_Atom::SetRadius(double radius)
{
    radius_ = radius;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void MD_Atom::Print(ostream &out)
{
}
