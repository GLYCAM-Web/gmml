#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include <vector>
#include <memory> // unique_ptr

#include "includes/Abstract/absResidue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

namespace cds
{
class cdsMolecule;
class cdsAtom;
class cdsResidue : public Abstract::absResidue, public glygraph::Node<cdsResidue>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    cdsResidue();
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() {return number_;}
    std::vector<cdsAtom*> getAtoms();
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
    std::vector<std::unique_ptr<cdsAtom>> atoms_;
    int number_;
};
}


#endif // RESIDUE_HPP
