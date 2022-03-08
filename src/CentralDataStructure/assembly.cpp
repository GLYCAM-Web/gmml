#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"

using cds::Assembly;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
Assembly::Assembly() : number_(0) {}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<cds::Molecule*> Assembly::getMolecules()
{
    std::vector<Molecule*> molecules;
    for(auto &molPtr : molecules_)
    {
        molecules.push_back(molPtr.get());
    }
    return molecules;
}

std::vector<cds::Residue*> Assembly::getResidues()
{
    std::vector<Residue*> residues;
    for(auto &molPtr : molecules_)
    {
        std::vector<Residue*> currentMoleculeResidues = molPtr->getResidues();
        residues.insert(residues.end(),
                std::make_move_iterator(currentMoleculeResidues.begin()),
                std::make_move_iterator(currentMoleculeResidues.end()) );
    }
    return residues;
}

std::vector<cds::Atom*> Assembly::getAtoms()
{
    std::vector<Atom*> atoms;
    for(auto &residue : this->getResidues())
    {
        std::vector<Atom*> currentResidueAtoms = residue->getAtoms();
        atoms.insert( atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but that's ok here.
                std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
