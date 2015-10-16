#ifndef ONTOLOGYVOCABULARY_HPP
#define ONTOLOGYVOCABULARY_HPP

#include <string>

namespace Ontology
{
    const std::string ONT_DOMAIN = "<http://gmmo.uga.edu#";
    const std::string ENTITY_COMMENT = "###  http://gmmo.uga.edu#";

    const std::string TYPE = "rdf:type";
    const std::string NAMED_INDIVIDUAL = "owl:NamedIndividual";

    /* Classes */
    const std::string Atom = ":Atom";
    const std::string RingAtom = ":RingAtom";
    const std::string SideAtom = ":SideAtom";
    const std::string Linkage = ":Linkage";
    const std::string Monosaccharide = ":Monosaccharide";
    const std::string Oligosaccharide = ":Oligosaccharide";
    const std::string PDB = ":PDB";
    const std::string Residue = ":Residue";
    const std::string SugarName = ":SugarName";

    /* Object Properties */
    const std::string hasAtom = ":hasAtom";
    const std::string hasChild = ":hasChild";
    const std::string hasChildAtomIndex = ":hasChildAtomIndex";
    const std::string hasGlycosidicLinkage = ":hasGlycosidicLinkage";
    const std::string hasNeighbor = ":hasNeighbor";
    const std::string hasOligo = ":hasOligo";
    const std::string hasParent = ":hasParent";
    const std::string hasParentAtomLinkage = ":hasParentAtomLinkage";
    const std::string hasResidue = ":hasResidue";
    const std::string hasRingAtom = ":hasRingAtom";
    const std::string hasRoot = ":hasRoot";
    const std::string hasSideAtom = ":hasSideAtom";
    const std::string hasSugarName = ":hasSugarName";
    const std::string hasmemberOfResidue = ":memberOfResidue"; ///remove maybe?

    /* Datatype Properties */
    const std::string anomeric_status = ":anomeric_status";
    const std::string chemical_code_str = ":chemical_code_str";
    const std::string configuration = ":configuration";
    const std::string cycle_atom_str = ":cycle_atom_str";
    const std::string derivative = ":derivative";
    const std::string glycosidic_linkage_str = ":glycosidic_linkage_str";
    const std::string group_name = ":group_name"; ///remove?
    const std::string id = ":id";
    const std::string isomer = ":isomer";
    const std::string linkage_str = ":linkage_str";
    const std::string mono_name = ":monosaccharide_name";
    const std::string mono_short_name = ":monosaccharide_short_name";
    const std::string mono_stereo_name = ":monosaccharide_stereochem_name";
    const std::string mono_stereo_short_name = ":monosaccharide_stereochem_short_name";
    const std::string oligo_name = ":oligo_name";
    const std::string orientation = ":orientation";
    const std::string path = ":path";
    const std::string ring_index = ":ring_index";
    const std::string ring_type = ":ring_type";
    const std::string side_index = ":side_index";
    const std::string x = ":x_crd";
    const std::string y = ":y_crd";
    const std::string z = ":z_crd";
}

#endif // ONTOLOGYVOCABULARY_HPP
