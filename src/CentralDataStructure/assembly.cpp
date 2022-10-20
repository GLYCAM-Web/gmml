#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/numbers.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp" // bondIfClose
#include "includes/CentralDataStructure/cdsFunctions.hpp"

using cds::Assembly;
using cds::Molecule;
using cds::Residue;
using cds::Atom;

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<Molecule*> Assembly::getMolecules() const
{
    std::vector<Molecule*> molecules;
    for(auto &molPtr : molecules_)
    {
        molecules.push_back(molPtr.get());
    }
    return molecules;
}

std::vector<Residue*> Assembly::getResidues() const
{
    std::vector<Residue*> residues;
    for(auto &molPtr : this->getMolecules())
    {
        std::vector<Residue*> currentMoleculeResidues = molPtr->getResidues();
        residues.insert(residues.end(),
                std::make_move_iterator(currentMoleculeResidues.begin()),
                std::make_move_iterator(currentMoleculeResidues.end()) );
    }
    return residues;
}

std::vector<Atom*> Assembly::getAtoms() const
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
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void Assembly::addMolecule(std::unique_ptr<Molecule> myMolecule)
{
    molecules_.push_back(std::move(myMolecule));
}

const Atom* Assembly::findAtom(const int& serialNumber) const
{
    return codeUtils::findElementWithNumber(this->getAtoms(), serialNumber);
}
