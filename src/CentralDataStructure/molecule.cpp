#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"

using cds::Molecule;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
Molecule::Molecule() : number_(0) {}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////

std::vector<cds::Residue*> Molecule::getResidues()
{
    std::vector<Residue*> residues;
    for(auto &residuePtr : residues_)
    {
        residues.push_back(residuePtr.get());
    }
    return residues;
}

std::vector<cds::Atom*> Molecule::getAtoms()
{
    std::vector<Atom*> atoms;
    for(auto &residuePtr : residues_)
    {
        std::vector<Atom*> currentResidueAtoms = residuePtr->getAtoms();
        atoms.insert( atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but that's ok here.
                std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
