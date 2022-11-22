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
    Residue() {} //{std::cout << "Residue default ctor\n";}
    Residue(const std::string& residueName, const Residue *referenceResidue);
    virtual ~Residue() {} //std::cout << "cdsResidue Dtor for " << this->getName() << "\n";}
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() const {return number_;}
    inline virtual const std::string& getName() const {return name_;}
    std::vector<Atom*> getAtoms() const;
    std::vector<std::string> getAtomNames() const;
    std::string getId(std::string moleculeNumber = constants::sNotSet) const;
    std::vector<Coordinate*> getCoordinates() const;
    const Coordinate* getGeometricCenter();
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    inline void setName(const std::string& s) {name_ = s;}
    inline void setAtoms(std::vector<std::unique_ptr<Atom>> v) {atoms_ = std::move(v);}
    void addAtom(std::unique_ptr<Atom> myAtom);
    //void addAtom(Atom* myAtom);
    bool deleteAtom(const Atom* atom);
    std::vector<std::unique_ptr<Atom>> extractAtoms() {return std::move(atoms_);}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    typename std::vector<std::unique_ptr<Atom>>::iterator FindPositionOfAtom(const Atom* queryAtom);
    Atom* FindAtom(const std::string queryName) const;
    Atom* FindAtom(const int& queryNumber) const;
    bool contains(const Atom* queryAtom) const;
    std::vector<const Atom*> getAtomsConnectedToOtherResidues() const;
    void MakeDeoxy(std::string oxygenNumber);
    const Coordinate* calculateGeometricCenter();
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
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
    Coordinate geometricCenter_;
};
} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
