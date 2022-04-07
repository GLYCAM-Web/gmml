#include "includes/CentralDataStructure/cdsAtom.hpp"
#include "includes/CentralDataStructure/cdsResidue.hpp"
#include "includes/CentralDataStructure/cdsMolecule.hpp"

using cds::cdsMolecule;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
cdsMolecule::cdsMolecule() : number_(0) {}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<cds::cdsResidue*> cdsMolecule::getResidues()
{
    std::vector<cdsResidue*> residues;
    for(auto &residuePtr : residues_)
    {
        residues.push_back(residuePtr.get());
    }
    return residues;
}

std::vector<cds::cdsAtom*> cdsMolecule::getAtoms()
{
    std::vector<cdsAtom*> atoms;
    for(auto &residuePtr : residues_)
    { // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but it's scoped to here.
        std::vector<cdsAtom*> currentResidueAtoms = residuePtr->getAtoms();
        atoms.insert( atoms.end(),
                std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
