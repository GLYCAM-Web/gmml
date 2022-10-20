#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSMOLECULE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSMOLECULE_HPP

#include <vector>
#include <memory> // unique_ptr

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
#include "includes/CentralDataStructure/Writers/cdsOffWriter.hpp"

namespace cds
{
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
    std::vector<residueT*> getResidues() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int i) {number_ = i;}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    void addResidue(std::unique_ptr<residueT> myResidue);
    residueT* addResidue(residueT* myResidue);
    residueT* createNewResidue(const std::string& residueName, const residueT& positionReferenceResidue);
    typename std::vector<std::unique_ptr<residueT>>::iterator findPositionOfResidue(const residueT* queryResidue);
    std::vector<residueT*> getResidues(std::vector<std::string> queryNames);
    residueT* getResidue(const std::string& queryName);
    void deleteResidue(residueT*);
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
    void WritePdb(std::ostream& stream) const;
    void WriteOff(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<residueT>> residues_;
    int number_;
};
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
//template< class residueT, class atomT >
//std::vector<const residueT*> cdsMolecule<residueT, atomT>::getResidues() const
//{
//    std::vector<const residueT*> residues;
//    for(auto &residuePtr : residues_)
//    {
//        residues.push_back(residuePtr.get());
//    }
//    return residues;
//}

template< class residueT, class atomT >
std::vector<residueT*> cdsMolecule<residueT, atomT>::getResidues() const
{
    std::vector<residueT*> residues;
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
//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
//template< class residueT, class atomT>
//void cdsMolecule<residueT, atomT>::addResidue(const residueT& residue)
//{ // This is bad: residue contains a vector of unique_ptr, so you don't want to copy that.
//    residues_.push_back(std::make_unique<residueT>(residue));
//}

template< class residueT, class atomT>
void cdsMolecule<residueT, atomT>::addResidue(std::unique_ptr<residueT> myResidue)
{ // This is good: myResidue contains a vector of unique_ptr, so you don't want to copy that.
    residues_.push_back(std::move(myResidue));
}

template< class residueT, class atomT>
residueT* cdsMolecule<residueT, atomT>::addResidue(residueT* myResidue)
{
    residues_.push_back(std::make_unique<residueT>(myResidue));
    return residues_.back();
}

template< class residueT, class atomT >
residueT* cdsMolecule<residueT, atomT>::createNewResidue(const std::string& residueName, const residueT& positionReferenceResidue)
{
    //Where the residue is in the vector matters. It should go after the reference residue.
    auto position = this->findPositionOfResidue(&positionReferenceResidue);
    if (position != residues_.end())
    {
        ++position; // it is ok to insert at end(). I checked. It was ok. Ok.
        position = residues_.insert(position, std::make_unique<residueT>(residueName, &positionReferenceResidue));
        gmml::log(__LINE__,__FILE__,gmml::INF, "New residue named " + residueName + " has been born; You're welcome.");
    }
    else
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create residue named " + residueName + " as referenceResidue was not found\n");
    }
    return (*position).get(); // Wow ok, so dereference the reference to a uniquePtr, then use get() to create a raw ptr.
}

template< class residueT, class atomT >
typename std::vector<std::unique_ptr<residueT>>::iterator cdsMolecule<residueT, atomT>::findPositionOfResidue(const residueT* queryResidue)
{
    auto i = residues_.begin();
    auto e = residues_.end();
    while (i != e)
    {
        if (queryResidue == i->get())
        {
            return i;
        }
        else
        {
            ++i;
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find position of " + queryResidue->getName() + " in vector\n"); // every class should have a print?
    return e;
}

template< class residueT, class atomT >
typename std::vector<residueT*> cdsMolecule<residueT, atomT>::getResidues(std::vector<std::string> queryNames)
{
	return codeUtils::getElementsWithNames(this->getResidues(), queryNames);
}

template< class residueT, class atomT >
residueT* cdsMolecule<residueT, atomT>::getResidue(const std::string& queryName)
{
    return codeUtils::findElementWithName(this->getResidues(), queryName);
}

template< class residueT, class atomT >
void cdsMolecule<residueT, atomT>::deleteResidue(residueT* residue)
{
    std::cout << "Gonna erase this mofo: " << std::endl;
    auto i = this->findPositionOfResidue(residue); // auto makes my life easier
    if (i != residues_.end())
    {
        gmml::log(__LINE__,__FILE__,gmml::INF, "Residue " + residue->getName() + " has been erased. You're welcome.");
        i = residues_.erase(i);
    }
    std::cout << "Done " << std::endl;
    return;
}

//////////////////////////////////////////////////////////
//                    DISPLAY                           //
//////////////////////////////////////////////////////////
template< class residueT, class atomT >
void cdsMolecule<residueT, atomT>::WritePdb(std::ostream& stream) const
{
    for (auto &residue : this->getResidues())
    {
        residue->WritePdb(stream);
    }
    stream << "TER\n";
}

template< class residueT, class atomT >
void cdsMolecule<residueT, atomT>::WriteOff(std::ostream& stream) const
{
	cds::WriteMoleculeToOffFile(this->getResidues(), stream, this->getName());
}

} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP
