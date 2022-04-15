#include "../../includes/MolecularMetadata/element.hpp"

using MolecularMetadata::Element;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Element::Element() {}
Element::~Element(){}
Element ::Element(Element *element)
{
    this->symbol_=element->GetSymbol();
    this->name_=element->GetName();
    this->mass_=element->GetMass();
    this->exact_mass_=element->GetExactMass();
    this->ionization_energy_=element->GetIonizationEnergy();
    this->election_affinity_=element->GetElectionAffinity();
    this->election_negativity_=element->GetElectionNegativity();
    this->covalent_radius_=element->GetCovalentRadius();
    this->van_der_waals_radius_=element->GetVanDerWaalsRadius();
    this->boiling_point_=element->GetBoilingPoint();
    this->melting_point_=element->GetMeltingPoint();
}
Element :: Element(Element& element)
{
    this->symbol_=element.GetSymbol();
    this->name_=element.GetName();
    this->mass_=element.GetMass();
    this->exact_mass_=element.GetExactMass();
    this->ionization_energy_=element.GetIonizationEnergy();
    this->election_affinity_=element.GetElectionAffinity();
    this->election_negativity_=element.GetElectionNegativity();
    this->covalent_radius_=element.GetCovalentRadius();
    this->van_der_waals_radius_=element.GetVanDerWaalsRadius();
    this->boiling_point_=element.GetBoilingPoint();
    this->melting_point_=element.GetMeltingPoint();
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string Element::GetSymbol()
{
    return symbol_;
}
std::string Element::GetName()
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
    return election_negativity_;
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
void Element::SetSymbol(std::string symbol)
{
    symbol_ = symbol;
}

void Element::SetName(std::string name)
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
void Element::Print(std::ostream &out)
{

    out << "*******Printing Element :" << name_ << " *******" << std::endl;
    out << "Element symbol: " << symbol_<< std::endl;
    out << "Element name: " << name_<< std::endl;
    out << "Element mass: " << mass_<< std::endl;
    out << "Element exact mass: " << exact_mass_<< std::endl;
    out << "Element ionization energy: " << ionization_energy_<< std::endl;
    out << "Element election affinity: " << election_affinity_<< std::endl;
    out << "Element election negativity: " << election_negativity_<< std::endl;
    out << "Element covalent radius: " << covalent_radius_<< std::endl;
    out << "Element van der waals radius: " << van_der_waals_radius_<< std::endl;
    out << "Element boiling point: " << boiling_point_<< std::endl;
    out << "Element melting point: " << melting_point_<< std::endl;
}
