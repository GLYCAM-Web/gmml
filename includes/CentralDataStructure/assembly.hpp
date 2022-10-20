#ifndef INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include <string>
#include <vector>
#include <memory> // unique_ptr
#include <algorithm> // std::find

namespace cds
{
class Assembly : public glygraph::Node<Assembly>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    Assembly() : number_(0) {}
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber()const {return number_ ;}
//    std::vector<const Atom*> getAtoms() const;
//    std::vector<const Residue*> getResidues() const;
//    std::vector<const Molecule*> getMolecules() const;
    std::vector<Atom*> getAtoms() const;
    std::vector<Residue*> getResidues() const;
    std::vector<Molecule*> getMolecules() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    void addMolecule(const Molecule& molecule);
    void addMolecule(std::unique_ptr<Molecule> myMolecule);
    const Atom* findAtom(const int& serialNumber) const;
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Molecule>> molecules_;
    int number_;
};
} // namespace
#endif // ASSEMBLY_HPP
