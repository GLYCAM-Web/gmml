#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_BEADRESIDUES_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_BEADRESIDUES_HPP

// Just planning to have these temporarily in the new central data structure while converting GP builder over.
// Think they can be replaced as I want to use distance instead of overlap, and can use the center of a residue anyway.
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
using cds::Assembly;
using cds::Residue;
using cds::Atom;
namespace beads
{
void Remove_Beads(Assembly &glycoprotein);
std::vector<Atom*> Add_Beads_To_Glycan(std::vector<Residue*> glycan_residues);
std::vector<Atom*> Add_Beads_To_Protein(Assembly &assembly);
}
#endif // GMML_INCLUDES_INTERNALPROGRAMS_BEADRESIDUES_HPP
