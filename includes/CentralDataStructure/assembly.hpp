#ifndef INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP

#include <string>
#include <vector>
#include <memory> // unique_ptr

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

namespace cds
{
class Atom;
class Residue;
class Molecule;
class Assembly : public glygraph::Node<Assembly>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    Assembly();
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() {return number_ ;}
    std::vector<Atom*> getAtoms();
    std::vector<Residue*> getResidues();
    std::vector<Molecule*> getMolecules();
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int& i) {number_ = i;}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
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
}

#endif // ASSEMBLY_HPP
