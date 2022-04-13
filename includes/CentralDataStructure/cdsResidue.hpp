#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include <vector>
#include <memory> // unique_ptr

#include "includes/Abstract/absResidue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/CodeUtils/logging.hpp"

namespace cds
{
class cdsAtom;
class cdsCoordinate;
template<class T>
class cdsResidue : public Abstract::absResidue, public glygraph::Node<cdsResidue<T>>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    cdsResidue() : number_(0) {}
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() {return number_;}
    std::vector<const T*> getAtoms() const;
    std::vector<std::string> getAtomNames() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    void CreateAtom(const std::string atomName, cdsCoordinate& atomCoord);
    void CreateAtom(T atom);
    bool DeleteAtom(T* atom);
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    typename std::vector<std::unique_ptr<T>>::iterator FindPositionOfAtom(const T* queryAtom) const;
    T* FindAtom(const std::string& queryName) const;
    T* FindAtom(const int& serialNumber) const;
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<T>> atoms_;
    int number_;
};

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
template< class T >
std::vector<const T*> cdsResidue<T>::getAtoms() const
{
    std::vector<const T*> atoms;
    for(auto &atomPtr : atoms_)
    {
        atoms.push_back(atomPtr.get());
    }
    return atoms;
}

template< class T >
std::vector<std::string> cdsResidue<T>::getAtomNames() const
{
    std::vector<std::string> foundAtomNames;
    for(auto &atom : this->getAtoms())
    {
        foundAtomNames.push_back(atom->getName());
    }
    return foundAtomNames;
}


//template< class T >
//const std::vector<const T*> cdsResidue<T>::getAtoms() const
//{
//    const std::vector<const T*> atoms;
//    for(auto &atomPtr : atoms_)
//    {
//        atoms.push_back(atomPtr.get());
//    }
//    return atoms;
//}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
template< class T >
void cdsResidue<T>::CreateAtom(const std::string atomName, cdsCoordinate& atomCoord)
{
    atoms_.push_back(std::make_unique<T>(atomName, this->GetName(), this->GetNumber(), this->GetInsertionCode(), atomCoord, this->GetChainId(), this->GetModelNumber()));
    return;
}

template< class T >
void cdsResidue<T>::CreateAtom(T atom)
{
    atoms_.push_back(std::make_unique<T>(atom));
    return;
}

template< class T >
bool cdsResidue<T>::DeleteAtom(T* atom)
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
template< class T >
typename std::vector<std::unique_ptr<T>>::iterator cdsResidue<T>::FindPositionOfAtom(const T* queryAtom) const
{
    auto i = atoms_.begin();
    auto e = atoms_.end();
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

template< class T >
T* cdsResidue<T>::FindAtom(const std::string& queryName) const
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

template< class T >
T* cdsResidue<T>::FindAtom(const int& serialNumber) const
{
    T* nullRecord = nullptr;
    for(auto &atomRecord : atoms_)
    {
        if (atomRecord->GetNumber() == serialNumber )
        {
            return atomRecord.get();
        }
    }
    return nullRecord;
}


} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
