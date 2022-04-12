#ifndef INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP

#include <string>
#include <vector>
#include <memory> // unique_ptr
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

namespace cds
{
class cdsAtom;
class cdsResidue;
class cdsMolecule;
template <class moleculeT, class residueT, class atomT>
class cdsAssembly : public glygraph::Node<cdsAssembly<moleculeT, residueT, atomT>>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    cdsAssembly() : number_(0) {}
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() {return number_ ;}
    std::vector<const atomT*> getAtoms() const;
    std::vector<const residueT*> getResidues() const;
    std::vector<const moleculeT*> getMolecules() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<moleculeT>> molecules_;
    int number_;
};
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
template <class moleculeT, class residueT, class atomT>
std::vector<const moleculeT*> cdsAssembly<moleculeT, residueT, atomT>::getMolecules() const
{
    std::vector<const moleculeT*> molecules;
    for(auto &molPtr : molecules_)
    {
        molecules.push_back(molPtr.get());
    }
    return molecules;
}

template <class moleculeT, class residueT, class atomT>
std::vector<const residueT*> cdsAssembly<moleculeT, residueT, atomT>::getResidues() const
{
    std::vector<const cdsResidue*> residues;
    for(auto &molPtr : molecules_)
    {
        std::vector<const cdsResidue*> currentMoleculeResidues = molPtr->getResidues();
        residues.insert(residues.end(),
                std::make_move_iterator(currentMoleculeResidues.begin()),
                std::make_move_iterator(currentMoleculeResidues.end()) );
    }
    return residues;
}

template <class moleculeT, class residueT, class atomT>
std::vector<const atomT*> cdsAssembly<moleculeT, residueT, atomT>::getAtoms() const
{
    std::vector<const cdsAtom*> atoms;
    for(auto &residue : this->getResidues())
    {
        std::vector<const cdsAtom*> currentResidueAtoms = residue->getAtoms();
        atoms.insert( atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but that's ok here.
                std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
} // namespace
#endif // ASSEMBLY_HPP
