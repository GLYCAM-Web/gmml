#ifndef INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ASSEMBLY_HPP

#include <string>
#include <vector>
#include <memory> // unique_ptr

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

namespace cds
{
class cdsAtom;
class cdsResidue;
class cdsMolecule;
class cdsAssembly : public glygraph::Node<cdsAssembly>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    cdsAssembly();
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() {return number_ ;}
    std::vector<cdsAtom*> getAtoms();
    std::vector<cdsResidue*> getResidues();
    std::vector<cdsMolecule*> getMolecules();
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
    std::vector<std::unique_ptr<cdsMolecule>> molecules_;
    int number_;
};
}

#endif // ASSEMBLY_HPP
