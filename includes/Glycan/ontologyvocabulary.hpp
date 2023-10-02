#ifndef ONTOLOGYVOCABULARY_HPP
#define ONTOLOGYVOCABULARY_HPP

#include <string>

namespace Ontology
{
    const std::string ONT_DOMAIN     = "<http://gmmo.uga.edu/#";
    const std::string ONT_PREFIX     = ":";
    const std::string ENTITY_COMMENT = " ###  http://gmmo.uga.edu#";

    const std::string TYPE  = "rdf:type";
    const std::string LABEL = "rdfs:label";

    const std::string TTL_FILE_PREFIX = "@prefix ttl: <http://www.semanticweb.org/owl/owlapi/turtle#> .\n"
                                        "@prefix owl: <http://www.w3.org/2002/07/owl#> .\n"
                                        "@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .\n"
                                        "@prefix xml: <http://www.w3.org/XML/1998/namespace> .\n"
                                        "@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .\n"
                                        "@prefix gmmo: <http://gmmo.uga.edu/#> .\n"
                                        "@prefix : <http://gmmo.uga.edu/#> .\n"
                                        "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .\n"
                                        "@base <http://gmmo.uga.edu/> .\n"
                                        "<http://gmmo.uga.edu/> rdf:type owl:Ontology .\n";

    // This is a silly way to organize this
    // todo - organize this by class?  Some things are shared?  Hmm

    /* Classes */
    const std::string Atom                 = ":Atom";
    const std::string RingAtom             = ":RingAtom";
    const std::string SideAtom             = ":SideAtom";
    const std::string Linkage              = ":Linkage";
    const std::string Monosaccharide       = ":Monosaccharide";
    const std::string Oligosaccharide      = ":Oligosaccharide";
    const std::string Terminal             = ":Terminal";
    const std::string Note                 = ":Note";
    const std::string PDB                  = ":PDB";
    const std::string PDBFile              = ":PDBFile";
    const std::string Residue              = ":Residue";
    const std::string SequenceResidue      = ":SequenceResidue";
    const std::string SugarName            = ":SugarName";
    const std::string ChemicalModification = ":ChemicalModification";

    /* Object Properties */
    const std::string hasAtom                  = ":hasAtom";
    const std::string hasChild                 = ":hasChild";
    const std::string hasChildMono             = ":hasChildMono";
    const std::string hasParentMono            = ":hasParentMono";
    const std::string hasChildAtomLinkage      = ":hasChildAtomLinkage";
    const std::string hasParentAtomLinkage     = ":hasParentAtomLinkage";
    const std::string hasGlycosidicLinkage     = ":hasGlycosidicLinkage";
    const std::string hasPhiAngle              = ":hasPhiAngle";
    const std::string hasPhiPrimeAngle         = ":hasPhiPrimeAngle";
    const std::string hasPsiAngle              = ":hasPsiAngle";
    const std::string hasOmegaAngle            = ":hasOmegaAngle";
    const std::string hasNeighbor              = ":hasNeighbor";
    const std::string hasOligo                 = ":hasOligo";
    const std::string hasMono                  = ":hasMono";
    const std::string hasTerminal              = ":hasTerminal";
    const std::string hasSequenceResidue       = ":hasSequenceResidue";
    const std::string hasResidueName           = ":hasResidueName";
    const std::string isConnectedTo            = ":isConnectedTo";
    const std::string hasNote                  = ":hasNote";
    const std::string hasParent                = ":hasParent";
    const std::string hasResidue               = ":hasResidue";
    const std::string hasRingAtom              = ":hasRingAtom";
    const std::string hasCore                  = ":hasCore";
    const std::string hasSideAtom              = ":hasSideAtom";
    const std::string hasSugarName             = ":hasSugarName";
    const std::string hasSNFGName              = ":hasSNFGName";
    const std::string hasIndex                 = ":hasIndex";
    const std::string hasNameIndex             = ":hasNameIndex";
    const std::string hasIUPACIndex            = ":hasIUPACIndex";
    const std::string hasOligoParent           = ":hasOligoParent";
    const std::string hasAuthorSNFGName        = ":hasAuthorSNFGName";
    const std::string hasTitle                 = ":hasTitle";
    const std::string hasAuthors               = ":hasAuthors";
    const std::string hasJournal               = ":hasJournal";
    const std::string hasDOI                   = ":hasDOI";
    const std::string hasPMID                  = ":hasPMID";
    const std::string hasResolution            = ":hasResolution";
    const std::string hasBFactor               = ":hasBFactor";
    const std::string hasProteinID             = ":hasProteinID";
    const std::string hasFormula               = ":hasFormula";
    const std::string isChemicallyModified     = ":isChemicallyModified";
    const std::string hasModifiedTerminal      = ":hasModifiedTerminal";
    const std::string isNucleotide             = ":isNucleotide";
    const std::string isSaccharide             = ":isSaccharide";
    const std::string isProtein                = ":isProtein";
    const std::string isAttachedToProtein      = ":isAttachedToProtein";
    const std::string isNGlycan                = ":isNGlycan";
    const std::string isOGlycan                = ":isOGlycan";
    const std::string isCGlycan                = ":isCGlycan";
    const std::string isSGlycan                = ":isSGlycan";
    const std::string isPGlycan                = ":isPGlycan";
    const std::string hasRingO                 = ":hasRingO";
    const std::string hasRingN                 = ":hasRingN";
    const std::string hasNonRingO              = ":hasNonRingO";
    const std::string hasNonRingN              = ":hasNonRingN";
    const std::string anomericNeighborElement  = ":anomericNElement"; // redundant?
    const std::string anomericProperlyAssigned = ":anomericProperlyAssigned";

    /* Datatype Properties */
    const std::string input_file_path               = ":filePath";
    const std::string anomeric_status               = ":anomericStatus";
    const std::string stereochemistry_chemical_code = ":chemicalCode";
    const std::string configuration                 = ":configuration";
    const std::string ring_atoms                    = ":ringAtoms";
    const std::string derivative                    = ":derivative";
    const std::string seq_derivative                = ":sequenceDerivative";
    const std::string glycosidic_linkage            = ":glycosidicLinkage";
    const std::string linkageType                   = ":linkageType";
    const std::string id                            = ":identifier";
    const std::string residue_linkage               = ":residueLinkage";
    const std::string isomer                        = ":isomer";
    const std::string mono_name                     = ":monoName";
    const std::string author_mono_name              = ":authorMonoName";
    const std::string mono_short_name               = ":shortName";
    const std::string author_mono_short_name        = ":authorShortName";
    const std::string mono_stereo_name              = ":stereochemName";
    const std::string author_mono_stereo_name       = ":authorStereochemName";
    const std::string mono_stereo_short_name        = ":stereochemShortName";
    const std::string author_mono_stereo_short_name = ":authorStereochemShortName";
    const std::string oligo_name                    = ":oligoName";
    const std::string oligo_IUPAC_name              = ":oligoIUPACname";
    const std::string author_IUPAC_name             = ":authorIUPACName";
    const std::string oligo_sequence_name           = ":oligoSequenceName";
    const std::string oligo_residue_linkages        = ":oligoResidueLinks";
    const std::string note_type                     = ":NoteType";
    const std::string note_category                 = ":NoteCategory";
    const std::string note_description              = ":description";
    const std::string orientation                   = ":orientation";
    const std::string path                          = ":path";
    const std::string ring_index                    = ":ringIndex";
    const std::string ring_type                     = ":ringType";
    const std::string side_index                    = ":sideIndex";
    const std::string x                             = ":xCoordinate";
    const std::string y                             = ":yCoordinate";
    const std::string z                             = ":zCoordinate";
    const std::string coordinate                    = ":coordinate";
    const std::string BFMP                          = ":BFMP";
    const std::string fullBFMP                      = ":fullBFMP";
    const std::string glycosylationType             = ":glycosylationType";
    const std::string glycosylationResidue          = ":glycosylationResidue";
    const std::string glycosylationPair             = ":glycosylationPair";
    const std::string totalCHIEnergy                = ":totalCHIEnergy";
    const std::string phiCHIEnergy                  = ":phiCHIEnergy";
    const std::string psiCHIEnergy                  = ":psiCHIEnergy";
    const std::string omegaCHIEnergy                = ":omegaCHIEnergy";
    const std::string phiCHIFunction                = ":phiCHIFunction";
    const std::string psiCHIFunction                = ":psiCHIFunction";
    const std::string omegaCHIFunction              = ":omegaCHIFunction";

    /* SPARQL Query */
    const std::string PREFIX = "\nPREFIX : <http://gmmo.uga.edu/#>\n"
                               "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n"
                               "PREFIX owl: <http://www.w3.org/2002/07/owl#>\n"
                               "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>\n"
                               "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>\n";

    const std::string SELECT_CLAUSE      = "SELECT";
    const std::string SELECT_DISTINCT    = "SELECT DISTINCT";
    const std::string WHERE_CLAUSE       = "WHERE \n{\n";
    const std::string END_WHERE_CLAUSE   = "}\n";
    const std::string CSV_OUTPUT_FORMAT  = " \'Accept: text/csv\' ";
    const std::string JSON_OUTPUT_FORMAT = " \'Accept: application/json\' ";
    const std::string XML_OUTPUT_FORMAT  = " \'Accept: application/sparql-results+xml' ";
    const std::string DATA_STORE_ADDRESS =
        "http://192.168.1.52:8890/sparql"; /* "http://128.192.62.244:8890/sparql", "http://192.168.1.52:8890/sparql" */
    const std::string DATA_STORE_ADDRESS_GF = "http://gw_virt:8890/sparql";
    const std::string CURL_PREFIX           = "curl -g -s -H ";
    const std::string QUERY_PREFIX          = " --data-urlencode query=\'";
    const std::string QUERY_POSTFIX         = "\'";

} // namespace Ontology

#endif // ONTOLOGYVOCABULARY_HPP
