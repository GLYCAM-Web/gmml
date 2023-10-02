#include "includes/CentralDataStructure/ensemble.hpp"

using cds::Assembly;
using cds::Atom;
using cds::Ensemble;
using cds::Molecule;
using cds::Residue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
// Move Ctor
Ensemble::Ensemble(Ensemble&& other) noexcept : glygraph::Node<cds::Ensemble>(other)
{
    assemblies_ = std::move(other.assemblies_);
}

// Copy Ctor
Ensemble::Ensemble(const Ensemble& other) : glygraph::Node<cds::Ensemble>(other)
{
    for (auto& assembly : other.assemblies_)
    {
        assemblies_.push_back(std::make_unique<Assembly>((*assembly.get())));
    }
}

// Move and Copy assignment operator
Ensemble& Ensemble::operator=(Ensemble other)
{
    glygraph::Node<cds::Ensemble>::operator=(other); // ToDo ok?.
    swap(*this, other);
    return *this;
}

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<Atom*> Ensemble::getAtoms()
{
    std::vector<Atom*> atoms;
    for (auto& residue : this->getResidues())
    {
        std::vector<Atom*> currentResidueAtoms = residue->getAtoms();
        atoms.insert(atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but
                                  // that's ok here; it goes out of scope.
                     std::make_move_iterator(currentResidueAtoms.begin()),
                     std::make_move_iterator(currentResidueAtoms.end()));
    }
    return atoms;
}

std::vector<Residue*> Ensemble::getResidues()
{
    std::vector<Residue*> residues;
    for (auto& molPtr : this->getMolecules())
    {
        std::vector<Residue*> currentMoleculeResidues = molPtr->getResidues();
        residues.insert(residues.end(), std::make_move_iterator(currentMoleculeResidues.begin()),
                        std::make_move_iterator(currentMoleculeResidues.end()));
    }
    return residues;
}

std::vector<Molecule*> Ensemble::getMolecules() const
{
    std::vector<Molecule*> molecules;
    for (auto& assPtr : this->getAssemblies())
    {
        std::vector<Molecule*> currentMoleculeResidues = assPtr->getMolecules();
        molecules.insert(molecules.end(), std::make_move_iterator(currentMoleculeResidues.begin()),
                         std::make_move_iterator(currentMoleculeResidues.end()));
    }
    return molecules;
}

std::vector<Assembly*> Ensemble::getAssemblies() const
{
    std::vector<Assembly*> assemblies;
    for (auto& assPtr : assemblies_)
    {
        assemblies.push_back(assPtr.get()); // raw ptr from unique_ptr
    }
    return assemblies;
}

//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
void Ensemble::addAssembly(std::unique_ptr<Assembly> myAssembly)
{
    assemblies_.push_back(std::move(myAssembly));
    return;
}
