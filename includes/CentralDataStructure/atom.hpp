#ifndef INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP

#include "includes/CentralDataStructure/cdsAtom.hpp"
// Playing around with a generic concrete class that can be used when we don't need a more
// specialized atom class like for example prepAtom or pdbAtom.
// To be used by ParsedResidue.
namespace cds
{
class Atom : public cdsAtom<Atom>
{
};
} // namespace
#endif // ATOM_HPP
