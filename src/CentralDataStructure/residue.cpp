#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"

using cds::Residue;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
Residue::Residue() : number_(0) {}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<cds::Atom*> Residue::getAtoms()
{
    std::vector<Atom*> atoms;
    for(auto &atomPtr : atoms_)
    {
        atoms.push_back(atomPtr.get());
    }
    return atoms;
}
