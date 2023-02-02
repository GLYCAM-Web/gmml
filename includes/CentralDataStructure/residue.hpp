#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/CodeUtils/constants.hpp" // iNotSet

#include <vector>
#include <memory> // unique_ptr

using cds::Coordinate;

namespace cds
{
enum ResidueType {Protein, Sugar, Aglycone, Derivative, Solvent, Deoxy, Undefined};
class Residue : public glygraph::Node<Residue>
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
    inline virtual const std::string& getName() const {return name_;}
    std::vector<Atom*> getAtoms() const;
    std::vector<std::string> getAtomNames() const;
    std::string getId(std::string moleculeNumber = constants::sNotSet) const;
    std::vector<Coordinate*> getCoordinates() const;
    const Coordinate* getGeometricCenter();
    inline ResidueType GetType() const {return type_;}
    inline unsigned int getNumber() const {return number_;}
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setName(const std::string& s) {name_ = s;}
    inline void setAtoms(std::vector<std::unique_ptr<Atom>> v) {atoms_ = std::move(v);}
    void addAtom(std::unique_ptr<Atom> myAtom);
    //void addAtom(Atom* myAtom);
    bool deleteAtom(const Atom* atom);
    std::vector<std::unique_ptr<Atom>> extractAtoms() {return std::move(atoms_);}
    inline void SetType(ResidueType type) {type_ = type;}
    inline void setNumber(unsigned int i) {number_ = i;}
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
    ResidueType determineType(const std::string &residueName);
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
    virtual void Print(std::ostream& out) const;
    //////////////////////////////////////////////////////////
    //                  OPERATOR OVERLOADING                //
    //////////////////////////////////////////////////////////
    virtual bool operator == (const Residue& rhs) const { return (this->getIndex() == rhs.getIndex());}
    virtual bool operator != (const Residue& rhs) const { return (this->getIndex() != rhs.getIndex());}
    virtual bool operator > (const Residue& rhs) const { return (this->getNumber() > rhs.getNumber());}
    virtual bool operator < (const Residue& rhs) const { return (this->getNumber() < rhs.getNumber());}
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Atom>> atoms_;
    std::string name_ = "   ";
    Coordinate geometricCenter_;
    ResidueType type_ = Undefined;  // enum Type. See enum above.
    unsigned int number_ = 1; //constants::iNotSet; ToDo: For prep residues a default 1 value is good. Is there a reason not to?
};
} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
