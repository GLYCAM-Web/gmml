#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/Abstract/absResidue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

#include <vector>
#include <memory> // unique_ptr

namespace cds
{
class Residue : public Abstract::absResidue, public glygraph::Node<Residue>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    Residue() {}
    Residue(const std::string& residueName, const Residue *referenceResidue);
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() const {return number_;}
    inline virtual const std::string& getName() const {return name_;}
    std::vector<const Atom*> getAtoms() const;
    std::vector<Atom*> getAtoms();
    std::vector<std::string> getAtomNames() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    inline void setName(const std::string& s) {name_ = s;}
    void createAtom(const std::string atomName, Coordinate& atomCoord);
    void addAtom(std::unique_ptr<Atom> myAtom);
    bool deleteAtom(const Atom* atom);
    std::vector<std::unique_ptr<Atom>> extractAtoms() {return std::move(atoms_);}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    typename std::vector<std::unique_ptr<Atom>>::iterator FindPositionOfAtom(const Atom* queryAtom);
    const Atom* FindAtom(const std::string queryName) const;
    const Atom* FindAtom(const int& queryNumber) const;
    std::vector<const Atom*> getAtomsConnectedToOtherResidues() const;
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
    void WritePdb(std::ostream& stream, bool addTerCard = false) const;
    virtual void Print(std::ostream& out) const;
    //////////////////////////////////////////////////////////
    //                  OPERATOR OVERLOADING                //
    //////////////////////////////////////////////////////////
    virtual bool operator== ( Residue& rhs)  { return (this->getIndex() == rhs.getIndex());}
    virtual bool operator!= ( Residue& rhs)  { return (this->getIndex() != rhs.getIndex());}
    virtual bool operator> ( Residue& rhs)  { return (this->getNumber() > rhs.getNumber());}
    virtual bool operator< ( Residue& rhs)  { return (this->getNumber() < rhs.getNumber());}
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Atom>> atoms_;
    int number_ = 1;
    std::string name_ = "   ";
};
} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
