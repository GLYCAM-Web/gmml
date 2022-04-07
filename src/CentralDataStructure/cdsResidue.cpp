#include "includes/CentralDataStructure/cdsResidue.hpp"
#include "includes/CentralDataStructure/cdsAtom.hpp"

using cds::cdsResidue;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
cdsResidue::cdsResidue() : number_(0) {}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<cds::cdsAtom*> cdsResidue::getAtoms()
{
    std::vector<cdsAtom*> atoms;
    for(auto &atomPtr : atoms_)
    {
        atoms.push_back(atomPtr.get());
    }
    return atoms;
}
