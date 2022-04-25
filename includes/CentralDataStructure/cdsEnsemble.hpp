#ifndef INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP

#include <string>
#include <vector>
#include <memory> // unique_ptr
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

namespace cds
{
template <class assemblyT, class moleculeT, class residueT, class atomT>
class cdsEnsemble //: public glygraph::Node<cdsEnsemble<assemblyT, moleculeT, residueT, atomT>>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    std::vector<const atomT*> getAtoms() const;
    std::vector<const residueT*> getResidues() const;
    std::vector<const moleculeT*> getMolecules() const;
    std::vector<const assemblyT*> getAssemblies() const;
    std::vector<assemblyT*> getAssemblies();
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    void addAssembly(std::unique_ptr<assemblyT> myAssembly);
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
    std::vector<std::unique_ptr<assemblyT>> assemblies_;

};
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////

template <class assemblyT, class moleculeT, class residueT, class atomT>
std::vector<const assemblyT*> cdsEnsemble<assemblyT, moleculeT, residueT, atomT>::getAssemblies() const
{
    std::vector<const assemblyT*> assemblies;
    for(auto &assPtr : assemblies_)
    {
        assemblies.push_back(assPtr.get()); // raw ptr from unique_ptr
    }
    return assemblies;
}

template <class assemblyT, class moleculeT, class residueT, class atomT>
std::vector<assemblyT*> cdsEnsemble<assemblyT, moleculeT, residueT, atomT>::getAssemblies()
{
    std::vector<assemblyT*> assemblies;
    for(auto &assPtr : assemblies_)
    {
        assemblies.push_back(assPtr.get()); // raw ptr from unique_ptr
    }
    return assemblies;
}

template <class assemblyT, class moleculeT, class residueT, class atomT>
std::vector<const moleculeT*> cdsEnsemble<assemblyT, moleculeT, residueT, atomT>::getMolecules() const
{
    std::vector<const moleculeT*> molecules;
    for(auto &assPtr : this->getAssemblies())
    {
        std::vector<const moleculeT*> currentAssMols;
        currentAssMols.push_back(assPtr.getMolecules());
        molecules.insert(molecules.end(),
                std::make_move_iterator(currentAssMols.begin()),
                std::make_move_iterator(currentAssMols.end()) );
    }
    return molecules;
}

template <class assemblyT, class moleculeT, class residueT, class atomT>
std::vector<const residueT*> cdsEnsemble<assemblyT, moleculeT, residueT, atomT>::getResidues() const
{
    std::vector<const residueT*> residues;
    for(auto &molPtr : this->getMolecules())
    {
        std::vector<const residueT*> currentMoleculeResidues = molPtr->getResidues();
        residues.insert(residues.end(),
                std::make_move_iterator(currentMoleculeResidues.begin()),
                std::make_move_iterator(currentMoleculeResidues.end()) );
    }
    return residues;
}

template <class assemblyT, class moleculeT, class residueT, class atomT>
std::vector<const atomT*> cdsEnsemble<assemblyT, moleculeT, residueT, atomT>::getAtoms() const
{
    std::vector<const atomT*> atoms;
    for(auto &residue : this->getResidues())
    {
        std::vector<const atomT*> currentResidueAtoms = residue->getAtoms();
        atoms.insert( atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but that's ok here.
                std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
//template <class assemblyT, class moleculeT, class residueT, class atomT>
//void cdsEnsemble<assemblyT, moleculeT, residueT, atomT>::addAssembly(const assemblyT& assembly)
//{
//    assemblies_.push_back(std::make_unique<assemblyT>(assembly));
//    return;
//}
template <class assemblyT, class moleculeT, class residueT, class atomT>
void cdsEnsemble<assemblyT, moleculeT, residueT, atomT>::addAssembly(std::unique_ptr<assemblyT> myAssembly)
{
    assemblies_.push_back(std::move(myAssembly));
    return;
}

} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP
