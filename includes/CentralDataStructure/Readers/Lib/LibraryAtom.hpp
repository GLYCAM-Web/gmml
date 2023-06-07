#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYATOM_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_LIB_LIBRARYATOM_HPP_

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CodeUtils/constants.hpp"

namespace lib
{
class LibraryAtom  : public cds::Atom
{
public:
    LibraryAtom (const std::string& line);
    int getResidueIndex();
    int getAtomIndex();
    int getAtomicNumber(); // This should be a local overload of the one in Abs::Atom.
private:
    //////////////////////////////////////////////////////////
    //                         ATTRIBUTES                   //
    //////////////////////////////////////////////////////////
    int residue_index_ = constants::iNotSet;    // Residue index that the atom belongs to; Set by the 4th column of the atom section of a library file
    int atom_index_ = constants::iNotSet;       // Index of the atom in the belonging residue; Set by the 6th column of the atom section of a library file
    int atomic_number_ = constants::iNotSet;    // Atomic number of the atom; Set by the 7th column of the atom section of a library file
};
} // namespace
#endif
