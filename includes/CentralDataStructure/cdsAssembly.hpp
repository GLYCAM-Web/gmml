#ifndef INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP

#include <string>
#include <vector>
#include <memory> // unique_ptr
#include <algorithm> // std::find
//#include <thread> // thread

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp" // bondIfClose


namespace cds
{
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
    inline const int& getNumber()const {return number_ ;}
    std::vector<const atomT*> getAtoms() const;
    std::vector<const residueT*> getResidues() const;
    std::vector<const moleculeT*> getMolecules() const;
    std::vector<atomT*> getAtoms();
    std::vector<residueT*> getResidues();
    std::vector<moleculeT*> getMolecules();
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    void addMolecule(const moleculeT& molecule);
    void addMolecule(std::unique_ptr<moleculeT> myMolecule);
    const atomT* findAtom(const int& serialNumber) const;
    void bondAtomsByDistance();
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
std::vector<moleculeT*> cdsAssembly<moleculeT, residueT, atomT>::getMolecules()
{
    std::vector<moleculeT*> molecules;
    for(auto &molPtr : molecules_)
    {
        molecules.push_back(molPtr.get());
    }
    return molecules;
}

template <class moleculeT, class residueT, class atomT>
std::vector<const residueT*> cdsAssembly<moleculeT, residueT, atomT>::getResidues() const
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

template <class moleculeT, class residueT, class atomT>
std::vector<residueT*> cdsAssembly<moleculeT, residueT, atomT>::getResidues()
{
    std::vector<residueT*> residues;
    for(auto &molPtr : this->getMolecules())
    {
        std::vector<residueT*> currentMoleculeResidues = molPtr->getResidues();
        residues.insert(residues.end(),
                std::make_move_iterator(currentMoleculeResidues.begin()),
                std::make_move_iterator(currentMoleculeResidues.end()) );
    }
    return residues;
}

template <class moleculeT, class residueT, class atomT>
std::vector<const atomT*> cdsAssembly<moleculeT, residueT, atomT>::getAtoms() const
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

template <class moleculeT, class residueT, class atomT>
std::vector<atomT*> cdsAssembly<moleculeT, residueT, atomT>::getAtoms()
{
    std::vector<atomT*> atoms;
    for(auto &residue : this->getResidues())
    {
        std::vector<atomT*> currentResidueAtoms = residue->getAtoms();
        atoms.insert( atoms.end(), // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but that's ok here.
                std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
//template <class moleculeT, class residueT, class atomT>
//void cdsAssembly<moleculeT, residueT, atomT>::addMolecule(const moleculeT& molecule)
//{
//    molecules_.push_back(std::make_unique<moleculeT>(molecule));
//}

template <class moleculeT, class residueT, class atomT>
void cdsAssembly<moleculeT, residueT, atomT>::addMolecule(std::unique_ptr<moleculeT> myMolecule)
{
    molecules_.push_back(std::move(myMolecule));
}

template <class moleculeT, class residueT, class atomT>
const atomT* cdsAssembly<moleculeT, residueT, atomT>::findAtom(const int& serialNumber) const
{
	return codeUtils::findElementWithNumber(this->getAtoms(), serialNumber);
}

// Simple for now, asked P to make a fast one.
template <class moleculeT, class residueT, class atomT>
void cdsAssembly<moleculeT, residueT, atomT>::bondAtomsByDistance()
{
//    std::vector<std::thread> threads;
//    int maxThreads = 11;
//    int i = 1;
    std::vector<atomT*> atoms = this->getAtoms();
    //typename std::vector<std::unique_ptr<atomT>>::iterator;
    for(typename std::vector<atomT*>::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        atomT* atom1 = *it1;
        for(typename std::vector<atomT*>::iterator it2 = std::next(it1); it2 != atoms.end(); ++it2)
        {
            atomT* atom2 = *it2;
            //gmml::log(__LINE__,__FILE__,gmml::INF, "Checking " + atom1->GetId() + " vs " + atom2->GetId());
            // Need to mutex lock the add node function in atom1? What about accessing atom2 as a node in that context? Probably also needs to lock.
            //std::thread t(&atomicBonds::bondAtomsIfClose, atom1, atom2);
            //t.join(); This makes it linear.
            atomicBonds::bondAtomsIfClose(atom1, atom2);
        }
    }
    return;
}

} // namespace
#endif // ASSEMBLY_HPP
