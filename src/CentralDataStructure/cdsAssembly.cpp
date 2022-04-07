#include "includes/CentralDataStructure/cdsAssembly.hpp"

#include "../../includes/CentralDataStructure/cdsResidue.hpp"
#include "includes/CentralDataStructure/cdsAtom.hpp"
#include "includes/CentralDataStructure/cdsMolecule.hpp"

using cds::cdsAssembly;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
cdsAssembly::cdsAssembly() : number_(0) {}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<cds::cdsMolecule*> cdsAssembly::getMolecules()
{
    std::vector<cdsMolecule*> molecules;
    for(auto &molPtr : molecules_)
    {
        molecules.push_back(molPtr.get());
    }
    return molecules;
}

std::vector<cds::cdsResidue*> cdsAssembly::getResidues()
{
    std::vector<cdsResidue*> residues;
    for(auto &molPtr : molecules_)
    {
        std::vector<cdsResidue*> currentMoleculeResidues = molPtr->getResidues();
        residues.insert(residues.end(),
                std::make_move_iterator(currentMoleculeResidues.begin()),
                std::make_move_iterator(currentMoleculeResidues.end()) );
    }
    return residues;
}

std::vector<cds::cdsAtom*> cdsAssembly::getAtoms()
{
    std::vector<cdsAtom*> atoms;
    for(auto &residue : this->getResidues())
    {
        std::vector<cdsAtom*> currentResidueAtoms = residue->getAtoms();
        atoms.insert( atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but that's ok here.
                std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
