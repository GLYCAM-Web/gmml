#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPPERMUTATOR_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPPERMUTATOR_HPP_

#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06LinkageCodes.hpp"

// Reads a prep file like 0GB and generates all possible permutations of decorated residues like 2GB, 3GB, YGB etc.
// Going from 0GB with a charge of 0.1940 to 6GB with a charge of 0.0 requires adjusting the charge on the O6 from -0.688 to -0.458 and deleting the H60 that has a charge of 0.4240.
// Remember: aglycones have -0.194, and every extra branch requires a -0.194 charge to accommodate the extra 0.1940 from the extra non-reducing terminal (e.g. 0GB) residue.
// When making a substituted (aka decorated) residue, you need to do:
// adjusted charge on oxygen = "old oxygen charge" + "removed Hydrogen Charge" - 0.194.
namespace prep
{

std::vector<std::string> selectAllAtomsThatCanBeSubstituted(std::vector<cds::Atom*> queryAtoms);
std::vector<std::vector<std::string>> getCombinations(const std::vector<std::string>& elements);
std::vector<PrepResidue*> generatePrepPermuations(PrepResidue* starterResidue);

} //namespace
#endif /* INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPPERMUTATOR_HPP_ */
