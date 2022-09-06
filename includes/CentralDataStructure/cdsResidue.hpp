#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include <vector>
#include <memory> // unique_ptr

#include "includes/Abstract/absResidue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"

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
    inline virtual const std::string& getName() const {return name_;}
    std::vector<const atomT*> getAtoms() const;
    std::vector<atomT*> getAtoms();
    std::vector<std::string> getAtomNames() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    inline void setName(const std::string& s) {name_ = s;}
    void createAtom(const std::string atomName, Coordinate& atomCoord);
    void addAtom(std::unique_ptr<atomT> myAtom);
    bool deleteAtom(const atomT* atom);
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    typename std::vector<std::unique_ptr<atomT>>::iterator FindPositionOfAtom(const atomT* queryAtom);
    const atomT* FindAtom(const std::string queryName) const;
    const atomT* FindAtom(const int& queryNumber) const;
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
    void WritePdb(std::ostream& stream, bool addTerCard = false) const;
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
std::vector<atomT*> cdsResidue<atomT>::getAtoms()
{
    std::vector<atomT*> atoms;
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
//const std::vector<atomT*> cdsResidue<atomT>::getAtoms()
//{
//    std::vector<atomT*> atoms;
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
bool cdsResidue<atomT>::deleteAtom(const atomT* atom)
{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
    auto i = this->FindPositionOfAtom(atom); // auto makes my life easier
    if (i != atoms_.end())
    {
       i = atoms_.erase(i);
       gmml::log(__LINE__,__FILE__,gmml::INF, "Atom " + atom->getName() + " has been erased. You're welcome.");
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
    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find " + queryAtom->getName() + " in atom records\n");
    return e;
}

template< class atomT >
const atomT* cdsResidue<atomT>::FindAtom(const std::string queryName) const
{
	return codeUtils::findElementWithName(this->getAtoms(), queryName);
}

template< class atomT >
const atomT* cdsResidue<atomT>::FindAtom(const int& queryNumber) const
{
	return codeUtils::findElementWithNumber(this->getAtoms(), queryNumber);
}
//////////////////////////////////////////////////////////
//                    DISPLAY                           //
//////////////////////////////////////////////////////////
template< class atomT >
void cdsResidue<atomT>::WritePdb(std::ostream& stream, const bool addTerCard) const
{
    for(auto &atom : this->getAtoms())
    {
        atom->WritePdb(stream);
    }
    if(addTerCard)
    {
        stream << "TER\n";
    }
}

template< class atomT >
void cdsResidue<atomT>::Print(std::ostream& out) const
{
	out << this->getName();
	for (auto &atom : this->getAtoms())
	{
		atom->Print(out);
	}
}

} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
