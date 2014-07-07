#include <iomanip>
#include "../../../includes/FileSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyatompair.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAtom::TopologyAtom() {}

TopologyAtom::TopologyAtom(int atom_index, string atom_name, string type, double atom_charge, int atomic_number, double atom_mass, ExcludedAtomNames excluded_atoms,
                           int number_of_excluded_atoms, double radii, double screen, string tree_chain_classification, string residue_name) :
    index_(atom_index), atom_name_(atom_name), type_(type), atom_charge_(atom_charge), atomic_number_(atomic_number), atom_mass_(atom_mass), number_of_excluded_atoms_(number_of_excluded_atoms),
    radii_(radii), screen_(screen), tree_chain_classification_(tree_chain_classification), residue_name_(residue_name)
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
string TopologyAtom::GetAtomName()
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
vector<string> TopologyAtom::GetExcludedAtoms()
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
string TopologyAtom::GetTreeChainClassification()
{
    return tree_chain_classification_;
}
int TopologyAtom::GetNumberOfExcludedAtoms()
{
    return number_of_excluded_atoms_;
}
string TopologyAtom::GetType()
{
    return type_;
}
int TopologyAtom::GetIndex()
{
    return index_;
}
string TopologyAtom::GetResidueName()
{
    return residue_name_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAtom::SetAtomName(const string atom_name)
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
void TopologyAtom::SetExcludedAtoms(vector<string> excluded_atoms)
{
    excluded_atoms_.clear();
    for(vector<string>::iterator it = excluded_atoms.begin(); it != excluded_atoms.end(); it++)
    {
        excluded_atoms_.push_back(*it);
    }
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
void TopologyAtom::SetType(string type)
{
    type_ = type;
}
void TopologyAtom::SetIndex(int index)
{
    index_ = index;
}
void TopologyAtom::SetResidueName(string residue_name)
{
    residue_name_ = residue_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAtom::Print(ostream &out)
{
    out << "Index: " << index_ << endl
           << "Residue Name: " << residue_name_ << endl
        << "Atom Name: " << atom_name_ << endl
        << "Charge: " << scientific << setprecision(8) << atom_charge_  << endl
        << "Atomic Number: " << atomic_number_ << endl
        << "Mass: " << scientific << setprecision(8) << atom_mass_ << endl
        << "Excluded Atoms: ";
    for(vector<string>::iterator it = excluded_atoms_.begin(); it != excluded_atoms_.end(); it++)
    {
        string atom_name = (*it);
        out << atom_name << "; ";
    }
    out << endl
        << "Radii: " << scientific << setprecision(8) << radii_ << endl
        << "Screen: " << scientific << setprecision(8) << screen_ << endl
        << "Tree Chain Classification: " << tree_chain_classification_ << endl
        << "Number of Excluded Atoms: " << number_of_excluded_atoms_ << endl
        << "Atom Type: " << type_ << endl << endl;
}


