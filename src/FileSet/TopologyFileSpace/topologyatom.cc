
#include "../../../includes/FileSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyatomtype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAtom::TopologyAtom() {}

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
string TopologyAtom::GetAtomType()
{
    return atom_type_;
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
char TopologyAtom::GetTreeChainClassification()
{
    return tree_chain_classification_;
}
TopologyAtom::TopologyLennardJonesMap TopologyAtom::GetLennardJonesA()
{
    return lennard_jones_a_;
}
TopologyAtom::TopologyLennardJonesMap TopologyAtom::GetLennardJonesB()
{
    return lennard_jones_b_;
}
int TopologyAtom::GetNumberOfExcludedAtoms()
{
    return number_of_excluded_atoms_;
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
void TopologyAtom::SetTreeChainClasification(char tree_chain_classification)
{
    tree_chain_classification_ = tree_chain_classification;
}
void TopologyAtom::SetNumberOfExcludedAtoms(int number_of_excluded_atoms)
{
    number_of_excluded_atoms_ = number_of_excluded_atoms;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAtom::Print(ostream &out)
{
}


