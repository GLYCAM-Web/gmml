#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include <vector>
#include <memory> // unique_ptr

#include "includes/Abstract/absResidue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/CodeUtils/logging.hpp"

using GeometryTopology::Coordinate;
namespace cds
{
template<class atomT>
class cdsResidue : public Abstract::absResidue, public glygraph::Node<cdsResidue<atomT>>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    cdsResidue() {}
    cdsResidue(const std::string& residueName, const cdsResidue *referenceResidue);
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() const {return number_;}
    inline const std::string& getName() const {return name_;}
    std::vector<const atomT*> getAtoms() const;
    std::vector<std::string> getAtomNames() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    inline void setName(const std::string& s) {name_ = s;}
    void createAtom(const std::string atomName, Coordinate& atomCoord);
    void addAtom(std::unique_ptr<atomT> myAtom);
    bool deleteAtom(atomT* atom);
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    typename std::vector<std::unique_ptr<atomT>>::iterator FindPositionOfAtom(const atomT* queryAtom);
    atomT* FindAtom(const std::string queryName) const;
    atomT* FindAtom(const int& queryNumber) const;
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<atomT>> atoms_;
    int number_ = 1;
    std::string name_ = "   ";
};

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
template< class atomT >
cdsResidue<atomT>::cdsResidue(const std::string& residueName, const cdsResidue *referenceResidue)
{
    this->setName(residueName);
    this->setNumber(referenceResidue->getNumber() + 1);
}

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
template< class atomT >
std::vector<const atomT*> cdsResidue<atomT>::getAtoms() const
{
    std::vector<const atomT*> atoms;
    for(auto &atomPtr : atoms_)
    {
        atoms.push_back(atomPtr.get());
    }
    return atoms;
}

template< class atomT >
std::vector<std::string> cdsResidue<atomT>::getAtomNames() const
{
    std::vector<std::string> foundAtomNames;
    for(auto &atom : this->getAtoms())
    {
        foundAtomNames.push_back(atom->getName());
    }
    return foundAtomNames;
}


//template< class atomT >
//const std::vector<const atomT*> cdsResidue<atomT>::getAtoms() const
//{
//    const std::vector<const atomT*> atoms;
//    for(auto &atomPtr : atoms_)
//    {
//        atoms.push_back(atomPtr.get());
//    }
//    return atoms;
//}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
template< class atomT >
void cdsResidue<atomT>::createAtom(const std::string atomName, Coordinate& atomCoord)
{
//    atoms_.push_back(std::make_unique<atomT>(atomName, this->getName(), this->getNumber(), this->GetInsertionCode(), atomCoord, this->GetChainId(), this->GetModelNumber()));
    atoms_.push_back(std::make_unique<atomT>(atomName, atomCoord));
    return;
}

template< class atomT >
void cdsResidue<atomT>::addAtom(std::unique_ptr<atomT> myAtom)
{
    atoms_.push_back(std::move(myAtom));
    return;
}


template< class atomT >
bool cdsResidue<atomT>::deleteAtom(atomT* atom)
{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
    auto i = this->FindPositionOfAtom(atom); // auto makes my life easier
    if (i != atoms_.end())
    {
       i = atoms_.erase(i);
       gmml::log(__LINE__,__FILE__,gmml::INF, "Atom " + atom->GetId() + " has been erased. You're welcome.");
       return true;
    }
    return false;
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
template< class atomT >
typename std::vector<std::unique_ptr<atomT>>::iterator cdsResidue<atomT>::FindPositionOfAtom(const atomT* queryAtom)
{
    typename std::vector<std::unique_ptr<atomT>>::iterator i = atoms_.begin();
    typename std::vector<std::unique_ptr<atomT>>::iterator e = atoms_.end();
    while (i != e)
    {
        if (queryAtom == i->get())
        {
            return i;
        }
        else
        {
            ++i;
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find " + queryAtom->GetId() + " in atom records\n");
    return e;
}

template< class atomT >
atomT* cdsResidue<atomT>::FindAtom(const std::string queryName) const
{
    for(auto &atom : atoms_)
    {
        if(atom->GetName() == queryName)
        {
            return atom.get();
        }
    }
    return nullptr;
}

template< class atomT >
atomT* cdsResidue<atomT>::FindAtom(const int& queryNumber) const
{
    atomT* nullRecord = nullptr;
    for(auto &atom : atoms_)
    {
        if (atom->getNumber() == queryNumber )
        {
            return atom.get();
        }
    }
    return nullRecord;
}


} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
