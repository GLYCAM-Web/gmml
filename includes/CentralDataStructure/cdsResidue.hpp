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
template<class Tatom>
class cdsResidue : public Abstract::absResidue, public glygraph::Node<cdsResidue<Tatom>>
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
    std::vector<const Tatom*> getAtoms() const;
    std::vector<std::string> getAtomNames() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    void createAtom(const std::string atomName, Coordinate& atomCoord);
    void createAtom(Tatom atom);
    bool deleteAtom(Tatom* atom);
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    typename std::vector<std::unique_ptr<Tatom>>::iterator FindPositionOfAtom(const T* queryAtom) const;
    Tatom* FindAtom(const std::string queryName) const;
    Tatom* FindAtom(const int& serialNumber) const;
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Tatom>> atoms_;
    int number_;
};

//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
template< class Tatom >
std::vector<const Tatom*> cdsResidue<Tatom>::getAtoms() const
{
    std::vector<const Tatom*> atoms;
    for(auto &atomPtr : atoms_)
    {
        atoms.push_back(atomPtr.get());
    }
    return atoms;
}

template< class Tatom >
std::vector<std::string> cdsResidue<Tatom>::getAtomNames() const
{
    std::vector<std::string> foundAtomNames;
    for(auto &atom : this->getAtoms())
    {
        foundAtomNames.push_back(atom->getName());
    }
    return foundAtomNames;
}


//template< class Tatom >
//const std::vector<const Tatom*> cdsResidue<Tatom>::getAtoms() const
//{
//    const std::vector<const Tatom*> atoms;
//    for(auto &atomPtr : atoms_)
//    {
//        atoms.push_back(atomPtr.get());
//    }
//    return atoms;
//}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
template< class Tatom >
void cdsResidue<Tatom>::createAtom(const std::string atomName, Coordinate& atomCoord)
{
    atoms_.push_back(std::make_unique<Tatom>(atomName, this->GetName(), this->GetNumber(), this->GetInsertionCode(), atomCoord, this->GetChainId(), this->GetModelNumber()));
    return;
}

template< class Tatom >
void cdsResidue<Tatom>::createAtom(Tatom atom)
{
    atoms_.push_back(std::make_unique<Tatom>(atom));
    return;
}

template< class Tatom >
bool cdsResidue<Tatom>::deleteAtom(Tatom* atom)
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
template< class Tatom >
typename std::vector<std::unique_ptr<Tatom>>::iterator cdsResidue<Tatom>::FindPositionOfAtom(const Tatom* queryAtom) const
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

template< class Tatom >
Tatom* cdsResidue<Tatom>::FindAtom(const std::string queryName) const
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

template< class Tatom >
Tatom* cdsResidue<Tatom>::FindAtom(const int& serialNumber) const
{
    Tatom* nullRecord = nullptr;
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
