#include <iomanip>
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../includes/common.hpp"

using TopologyFileSpace::TopologyAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAtom::TopologyAtom()
{
    index_ = gmml::iNotSet;
    atom_name_ = "";
    type_ = "";
    atom_charge_ = gmml::dNotSet;
    atomic_number_ = gmml::iNotSet;
    atom_mass_ = gmml::dNotSet;
    excluded_atoms_ = ExcludedAtomNames();
    number_of_excluded_atoms_ = gmml::iNotSet;
    radii_ = gmml::dNotSet;
    screen_ = gmml::dNotSet;
    tree_chain_classification_ = "";
    residue_name_ = "";
}

TopologyAtom::TopologyAtom( int atom_index, std::string atom_name, std::string type,
                            double atom_charge, int atomic_number, double atom_mass,
                            ExcludedAtomNames excluded_atoms, int number_of_excluded_atoms,
                            double radii, double screen, std::string tree_chain_classification,
                            std::string residue_name) :
    atom_name_(atom_name), atom_charge_(atom_charge), atomic_number_(atomic_number),
    atom_mass_(atom_mass), radii_(radii), screen_(screen),
    tree_chain_classification_(tree_chain_classification),
    number_of_excluded_atoms_(number_of_excluded_atoms),
    type_(type), index_(atom_index), residue_name_(residue_name)
{
    excluded_atoms_.clear();
    for(ExcludedAtomNames::iterator it = excluded_atoms.begin(); it != excluded_atoms.end(); it++)
    {
        excluded_atoms_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string TopologyAtom::GetAtomName()
{
    return atom_name_;
}
double TopologyAtom::GetAtomCharge()
{
    return atom_charge_;
}
int TopologyAtom::GetAtomicNumber()
{
    return atomic_number_;
}
double TopologyAtom::GetAtomMass()
{
    return atom_mass_;
}
std::vector<std::string> TopologyAtom::GetExcludedAtoms()
{
    return excluded_atoms_;
}
double TopologyAtom::GetRadii()
{
    return radii_;
}
double TopologyAtom::GetScreen()
{
    return screen_;
}
std::string TopologyAtom::GetTreeChainClassification()
{
    return tree_chain_classification_;
}
int TopologyAtom::GetNumberOfExcludedAtoms()
{
    return number_of_excluded_atoms_;
}
std::string TopologyAtom::GetType()
{
    return type_;
}
int TopologyAtom::GetIndex()
{
    return index_;
}
std::string TopologyAtom::GetResidueName()
{
    return residue_name_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAtom::SetAtomName(const std::string atom_name)
{
    atom_name_ = atom_name;
}
void TopologyAtom::SetAtomCharge(double atom_charge)
{
    atom_charge_ = atom_charge;
}
void TopologyAtom::SetAtomicNumber(int atomic_number)
{
    atomic_number_ = atomic_number;
}
void TopologyAtom::SetAtomMass(double atom_mass)
{
    atom_mass_ = atom_mass;
}
void TopologyAtom::SetExcludedAtoms(std::vector<std::string> excluded_atoms)
{
    excluded_atoms_.clear();
    for(std::vector<std::string>::iterator it = excluded_atoms.begin(); it != excluded_atoms.end(); it++)
    {
        excluded_atoms_.push_back(*it);
    }
}
void TopologyAtom::AddExcludedAtom(std::string excluded_atom)
{
        excluded_atoms_.push_back(excluded_atom);
}
void TopologyAtom::SetRadii(double radii)
{
    radii_ = radii;
}
void TopologyAtom::SetScreen(double screen)
{
    screen_ = screen;
}
void TopologyAtom::SetTreeChainClasification(std::string tree_chain_classification)
{
    tree_chain_classification_ = tree_chain_classification;
}
void TopologyAtom::SetNumberOfExcludedAtoms(int number_of_excluded_atoms)
{
    number_of_excluded_atoms_ = number_of_excluded_atoms;
}
void TopologyAtom::SetType(std::string type)
{
    type_ = type;
}
void TopologyAtom::SetIndex(int index)
{
    index_ = index;
}
void TopologyAtom::SetResidueName(std::string residue_name)
{
    residue_name_ = residue_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAtom::Print(std::ostream &out)
{
    out << "Index: " << index_ << std::endl
           << "Residue Name: " << residue_name_ << std::endl
        << "Atom Name: " << atom_name_ << std::endl
        << "Charge: " << std::scientific << std::setprecision(8) << atom_charge_  << std::endl
        << "Atomic Number: " << atomic_number_ << std::endl
        << "Mass: " << std::scientific << std::setprecision(8) << atom_mass_ << std::endl
        << "Excluded Atoms: ";
    for(std::vector<std::string>::iterator it = excluded_atoms_.begin(); it != excluded_atoms_.end(); it++)
    {
        std::string atom_name = (*it);
        out << atom_name << "; ";
    }
    out << std::endl
        << "Radii: " << std::scientific << std::setprecision(8) << radii_ << std::endl
        << "Screen: " << std::scientific << std::setprecision(8) << screen_ << std::endl
        << "Tree Chain Classification: " << tree_chain_classification_ << std::endl
        << "Number of Excluded Atoms: " << number_of_excluded_atoms_ << std::endl
        << "Atom Type: " << type_ << std::endl << std::endl;
}
