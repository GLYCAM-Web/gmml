#ifndef ONTOLOGYVOCABULARY_HPP
#define ONTOLOGYVOCABULARY_HPP

#include <string>

namespace Ontology
{
    const std::string ONT_DOMAIN = "<http://gmmo.uga.edu/#";
    const std::string ONT_PREFIX = "gmmo:";
    const std::string ENTITY_COMMENT = " ###  http://gmmo.uga.edu#";

    const std::string TYPE = "rdf:type";
    const std::string LABEL = "rdfs:label";

    const std::string TTL_FILE_PREFIX = "@prefix : <http://www.semanticweb.org/owl/owlapi/turtle#> .\n"
                                        "@prefix owl: <http://www.w3.org/2002/07/owl#> .\n"
                                        "@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .\n"
                                        "@prefix xml: <http://www.w3.org/XML/1998/namespace> .\n"
                                        "@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .\n"
                                        "@prefix gmmo: <http://gmmo.uga.edu/#> .\n"
                                        "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .\n"
                                        "@base <http://gmmo.uga.edu/> .\n"
                                        "<http://gmmo.uga.edu/> rdf:type owl:Ontology .\n";

    /* Classes */
    const std::string Atom = "gmmo:Atom";
    const std::string RingAtom = "gmmo:RingAtom";
    const std::string SideAtom = "gmmo:SideAtom";
    const std::string Linkage = "gmmo:Linkage";
    const std::string monosaccharide = "gmmo:Monosaccharide";
    const std::string Oligosaccharide = "gmmo:Oligosaccharide";
    const std::string Terminal = "gmmo:Terminal";
    const std::string Note = "gmmo:Note";
    const std::string PDB = "gmmo:PDB";
    const std::string Residue = "gmmo:Residue";
    const std::string SequenceResidue = "gmmo:SequenceResidue";
    const std::string SugarName = "gmmo:SugarName";

    /* Object Properties */
    const std::string hasAtom = "gmmo:hasAtom";
    const std::string hasChild = "gmmo:hasChild";
    const std::string hasChildMono = "gmmo:hasChildMono";
    const std::string hasParentMono = "gmmo:hasParentMono";
    const std::string hasChildAtomLinkage = "gmmo:hasChildAtomLinkage";
    const std::string hasParentAtomLinkage = "gmmo:hasParentAtomLinkage";
    const std::string hasGlycosidicLinkage = "gmmo:hasGlycosidicLinkage";
    const std::string hasGlycosidicPhiAngle = "gmmo:hasGlycosidicPhiAngle";
    const std::string hasGlycosidicPsiAngle = "gmmo:hasGlycosidicPsiAngle";
    const std::string hasGlycosidicOmegaAngle = "gmmo:hasGlycosidicOmegaAngle";
    const std::string hasNeighbor = "gmmo:hasNeighbor";
    const std::string hasOligo = "gmmo:hasOligo";
    const std::string hasMono = "gmmo:hasMono";
    const std::string hasTerminal = "gmmo:hasTerminal";
    const std::string hasSequenceResidue = "gmmo:hasSequenceResidue";
    const std::string isConnectedTo = "gmmo:isConnectedTo";
    const std::string hasNote = "gmmo:hasNote";
    const std::string hasParent = "gmmo:hasParent";
    const std::string hasResidue = "gmmo:hasResidue";
    const std::string hasRingAtom = "gmmo:hasRingAtom";
    const std::string hasCore = "gmmo:hasCore";
    const std::string hasSideAtom = "gmmo:hasSideAtom";
    const std::string hasSugarName = "gmmo:hasSugarName";
    const std::string hasSNFGName = "gmmo:hasSNFGName";
    const std::string hasIndex = "gmmo:hasIndex";
    const std::string hasNameIndex = "gmmo:hasNameIndex";
    const std::string hasOligoParent = "gmmo:hasOligoParent";
    const std::string hasAuthorSNFGName = "gmmo:hasAuthorSNFGName";
    const std::string hasTitle = "gmmo:hasTitle";
    const std::string hasAuthors = "gmmo:hasAuthors";
    const std::string hasJournal = "gmmo:hasJournal";
    const std::string hasDOI = "gmmo:hasDOI";
    const std::string hasPMID = "gmmo:hasPMID";
    const std::string hasResolution = "gmmo:hasResolution";
    const std::string hasBFactor = "gmmo:hasBFactor";
    const std::string hasProteinID = "gmmo:hasProteinID";
    const std::string hasFormula = "gmmo:hasFormula";
    const std::string isNucleotide = "gmmo:isNucleotide";

    /* Datatype Properties */
    const std::string input_file_path = "gmmo:filePath";
    const std::string anomeric_status = "gmmo:anomericStatus";
    const std::string stereochemistry_chemical_code = "gmmo:stereochemistryChemicalCode";
    const std::string configuration = "gmmo:configuration";
    const std::string ring_atoms = "gmmo:ringAtoms";
    const std::string derivative = "gmmo:derivative";
    const std::string seq_derivative = "gmmo:sequenceDerivative";
    const std::string glycosidic_linkage = "gmmo:glycosidicLinkage";
    const std::string linkageIndeces = "gmmo:linkageIndeces";
    const std::string id = "gmmo:identifier";
    const std::string sequence_linkage = "gmmo:sequenceLinkage";
    const std::string isomer = "gmmo:isomer";
    const std::string mono_name = "gmmo:monosaccharideName";
    const std::string author_mono_name = "gmmo:authorMonosaccharideName";
    const std::string mono_short_name = "gmmo:monosaccharideShortName";
    const std::string author_mono_short_name = "gmmo:authorMonosaccharideShortName";
    const std::string mono_stereo_name = "gmmo:monosaccharideStereochemName";
    const std::string author_mono_stereo_name = "gmmo:authorMonosaccharideStereochemName";
    const std::string mono_stereo_short_name = "gmmo:monosaccharideStereochemShortName";
    const std::string author_mono_stereo_short_name = "gmmo:authorMonosaccharideStereochemShortName";
    const std::string oligo_name = "gmmo:oligoName";
    const std::string oligo_IUPAC_name = "gmmo:oligoIUPACname";
    const std::string author_oligo_name = "gmmo:authorOligoName";
    const std::string oligo_sequence_name = "gmmo:oligoSequenceName";
    const std::string oligo_residue_linkages = "gmmo:oligoResidueLinks";
    const std::string oligo_b_factor = "gmmo:oligoBFactor";
    const std::string note_type = "gmmo:NoteType";
    const std::string note_category = "gmmo:NoteCategory";
    const std::string note_description = "gmmo:description";
    const std::string orientation = "gmmo:orientation";
    const std::string path = "gmmo:path";
    const std::string ring_index = "gmmo:ringIndex";
    const std::string ring_type = "gmmo:ringType";
    const std::string side_index = "gmmo:sideIndex";
    const std::string x = "gmmo:xCoordinate";
    const std::string y = "gmmo:yCoordinate";
    const std::string z = "gmmo:zCoordinate";
    const std::string coordinate = "gmmo:coordinate";
    const std::string bfmp_ring_conformation = "gmmo:BFMPRingConformation";

    /* SPARQL Query */
    const std::string PREFIX = "PREFIX : <http://gmmo.uga.edu/#>\nPREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\nPREFIX owl: <http://www.w3.org/2002/07/owl#>\nPREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>\nPREFIX xsd: <http://www.w3.org/2001/XMLSchema#>\n";
    const std::string SELECT_CLAUSE = "SELECT";
    const std::string WHERE_CLAUSE = "WHERE {\n";
    const std::string END_WHERE_CLAUSE = "}\n";
    const std::string CSV_OUTPUT_FORMAT = " \'Accept: text/csv\' ";
    const std::string JSON_OUTPUT_FORMAT = " \'Accept: application/json\' ";
    const std::string XML_OUTPUT_FORMAT = " \'Accept: application/sparql-results+xml' ";
    const std::string DATA_STORE_ADDRESS = "http://192.168.1.52:8890/sparql"; /* "http://128.192.62.244:8890/sparql", "http://192.168.1.52:8890/sparql" */
    const std::string DATA_STORE_ADDRESS_GF = "http://gw_virt:8890/sparql";
    const std::string CURL_PREFIX = "curl -g -s -H";
    const std::string QUERY_PREFIX = " --data-urlencode query=\'";
    const std::string QUERY_POSTFIX = "\'";
}

#endif // ONTOLOGYVOCABULARY_HPP
