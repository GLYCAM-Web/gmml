
#include "../../includes/MolecularModeling/element.hpp"

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Element::Element() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string Element::GetSymbol()
{
    return symbol_;
}
string Element::GetName()
{
    return name_;
}
double Element::GetMass()
{
    return mass_;
}
double Element::GetExactMass()
{
    return exact_mass_;
}
double Element::GetIonizationEnergy()
{
    return ionization_energy_;
}
double Element::GetElectionAffinity()
{
    return election_affinity_;
}
double Element::GetElectionNegativity()
{
    return election_affinity_;
}
double Element::GetCovalentRadius()
{
    return covalent_radius_;
}
double Element::GetVanDerWaalsRadius()
{
    return van_der_waals_radius_;
}
double Element::GetBoilingPoint()
{
    return boiling_point_;
}
double Element::GetMeltingPoint()
{
    return melting_point_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Element::SetSymbol(string symbol)
{
    symbol_ = symbol;
}

void Element::SetName(string name)
{
    name_ = name;
}
void Element::SetMass(double mass)
{
    mass_ = mass;
}
void Element::SetExactMass(double exact_mass)
{
    exact_mass_ = exact_mass;
}
void Element::SetIonizationEnergy(double ionization_energy)
{
    ionization_energy_ = ionization_energy;
}
void Element::SetElectionAffinity(double election_affinity)
{
    election_affinity_ = election_affinity;
}
void Element::SetElectionNegativity(double election_negativity)
{
    election_negativity_ = election_negativity;
}
void Element::SetCovalentRadius(double covalent_radius)
{
    covalent_radius_ = covalent_radius;
}
void Element::SetVanDerWaalsRadius(double van_der_waals_radius)
{
    van_der_waals_radius_ = van_der_waals_radius;
}
void Element::SetBoilingPoint(double boiling_point)
{
    boiling_point_ = boiling_point;
}
void Element::SetMeltingPoint(double melting_point)
{
    melting_point_ = melting_point;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Element::Print(ostream &out)
{
}

