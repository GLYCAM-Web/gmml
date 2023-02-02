#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/Abstract/absResidue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

#include <vector>
#include <memory> // unique_ptr
#include <algorithm> //swap

using cds::Coordinate;

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
    Residue(Residue&& other) noexcept;     // Move Ctor
    Residue(const Residue& other);         // Copy Ctor
    Residue& operator=(Residue other);     // Move & Copy assignment operator.
    virtual ~Residue() {}                  // Dtor, virtual so that derived classes dtors will get triggered if possible.
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
    virtual bool operator == (const Residue& rhs) const { return (this->getIndex() == rhs.getIndex());}
    virtual bool operator != (const Residue& rhs) const { return (this->getIndex() != rhs.getIndex());}
    virtual bool operator > (const Residue& rhs) const { return (this->getNumber() > rhs.getNumber());}
    virtual bool operator < (const Residue& rhs) const { return (this->getNumber() < rhs.getNumber());}
    //////////////////////////////////////////////////////////
    //               Copy-Swap                    //
    //////////////////////////////////////////////////////////
    friend void swap(Residue& lhs, Residue& rhs)
    {
        std::cout << "Swapping" << std::endl;
        using std::swap;
        swap(lhs.atoms_, rhs.atoms_);
        swap(lhs.number_, rhs.number_);
        swap(lhs.name_, rhs.name_);
        swap(lhs.geometricCenter_, rhs.geometricCenter_);
        std::cout << "Swapped" << std::endl;
    }
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
