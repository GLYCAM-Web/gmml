#ifndef GMML_INCLUDES_INTERNALPROGRAMS_BEADRESIDUES_HPP
#define GMML_INCLUDES_INTERNALPROGRAMS_BEADRESIDUES_HPP

#include "includes/MolecularModeling/assembly.hpp"
#include "includes/MolecularModeling/residue.hpp"
#include "includes/MolecularModeling/atom.hpp"

namespace beads
{
    void Remove_Beads(MolecularModeling::Assembly& glycoprotein);
    AtomVector Add_Beads_To_Glycan(ResidueVector glycan_residues);
    AtomVector Add_Beads_To_Protein(MolecularModeling::Assembly& assembly);
} // namespace beads
#endif // GMML_INCLUDES_INTERNALPROGRAMS_BEADRESIDUES_HPP
