#include "../../includes/MolecularModeling/Molecule.hpp"


using namespace std;
using namespace MolecularModeling;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Molecule::Molecule(){

    cout<<"Molecule class is called"<<endl;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int Molecule::GetMoleculeIndex()
{
    return molecule_index_;
}

Molecule::AtomVector Molecule::GetMoleculeAtoms()
{
    return molecule_atoms_;
}

Molecule::ResidueVector Molecule::GetMoleculeResidues()
{
    return molecule_residues_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

void Molecule::SetMoleculeIndex(int molecule_index)
{
    molecule_index_= molecule_index;
}

void Molecule::SetMoleculeAtoms(AtomVector molecule_atoms)
{
    molecule_atoms_.clear();
    for(AtomVector::iterator it = molecule_atoms.begin(); it != molecule_atoms.end(); it++)
        molecule_atoms_.push_back(*it);
}

void Molecule::SetMoleculeResidues(ResidueVector molecule_residues)
{
    molecule_residues_.clear();
    for(ResidueVector::iterator it = molecule_residues.begin(); it != molecule_residues.end(); it++)
        molecule_residues_.push_back(*it);
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Molecule::Print(ostream &out)
{
        out<<"Printing molecule details"<<endl;
}

