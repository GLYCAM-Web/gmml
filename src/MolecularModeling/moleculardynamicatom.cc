#include "../../includes/MolecularModeling/moleculardynamicatom.hpp"
#include "../../includes/common.hpp"

using MolecularModeling::MolecularDynamicAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
MolecularDynamicAtom::MolecularDynamicAtom() : atom_type_(""), charge_(gmml::dNotSet), mass_(gmml::dNotSet), radius_(gmml::dNotSet) {}

MolecularDynamicAtom::MolecularDynamicAtom(MolecularDynamicAtom& moleculardynamicatom){
    this->atom_type_=moleculardynamicatom.GetAtomType();
    this->charge_=moleculardynamicatom.GetCharge();
    this->mass_=moleculardynamicatom.GetMass();
    this->radius_=moleculardynamicatom.GetRadius();
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string MolecularDynamicAtom::GetAtomType()
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
void MolecularDynamicAtom::SetAtomType(std::string atom_type)
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
void MolecularDynamicAtom::Print(std::ostream &out)
{
     out << "------------------------ Atom Type :" << atom_type_ << " --------------------------" << std::endl;
     out << "Atom Charge :" << charge_ << std::endl;
     out << "Atom Mass :" << mass_ << std::endl;
     out << "Atom Radius :" << radius_ << std::endl;
}

void MolecularDynamicAtom::PrintHet(std::ostream &out)
{
    out << "";
}
