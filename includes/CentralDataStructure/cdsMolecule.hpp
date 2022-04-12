#ifndef INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP

#include <vector>
#include <memory> // unique_ptr
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

namespace cds
{
class cdsAtom;
class cdsResidue;
template <class residueT, class atomT>
class cdsMolecule : public glygraph::Node<cdsMolecule<residueT,atomT>>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    cdsMolecule() : number_(0) {}
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() {return number_;}
    std::vector<const atomT*> getAtoms() const;
    std::vector<const residueT*> getResidues() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int i) {number_ = i;}
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
    std::vector<std::unique_ptr<residueT>> residues_;
    int number_;
};
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
template< class residueT, class atomT >
std::vector<const residueT*> cdsMolecule<residueT, atomT>::getResidues() const
{
    std::vector<const residueT*> residues;
    for(auto &residuePtr : residues_)
    {
        residues.push_back(residuePtr.get());
    }
    return residues;
}

template< class residueT, class atomT >
std::vector<const atomT*> cdsMolecule<residueT, atomT>::getAtoms() const
{
    std::vector<const atomT*> atoms;
    for(auto &residuePtr : residues_)
    { // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but it's scoped to here.
        std::vector<const atomT*> currentResidueAtoms = residuePtr->getAtoms();
        atoms.insert( atoms.end(),
                std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP
