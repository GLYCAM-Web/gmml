#include "../../../../includes/MolecularModeling/assembly.hpp"

void MolecularModeling::Assembly::QueryOntology(std::string searchType, std::string searchTerm, std::string output_file_type)
{
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
    query << "?pdb ?oligo_sequence ?residue_links ?glycosidic_linkage ?title"
             "?resolution ?Mean_B_Factor ?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI"; 
    query << Ontology::WHERE_CLAUSE;
    query << "?oligo        :oligoName              ?oligo_sequence.\n";
    query << "?pdb_file     :identifier             ?pdb.\n";
    query << "OPTIONAL {";
    query << "?pdb_file     :hasTitle               ?title.\n";
    query << "?pdb_file     :hasResolution          ?resolution.\n";
    query << "?pdb_file     :hasAuthors             ?authors.\n";
    query << "?pdb_file     :hasJournal             ?journal.\n";
    query << "?pdb_file     :hasDOI                 ?DOI.\n";
    query << "?pdb_file     :hasPMID                ?PMID.\n";
    query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.\n";
    query << "?oligo        :oligoName              ?oligo_sequence.\n";
    query << "?oligo        :oligoBFactor           ?oligo_mean_B_Factor.\n";
    query << "?pdb_file     :hasOligo	              ?oligo.\n";
    query << "?oligo        :oligoResidueLinks      ?residue_links.\n";
    query << "?linkage      :hasParent 	            ?oligo.\n";
    query << "?linkage      :glycosidicLinkage      ?glycosidic_linkage.}\n";
    query << Ontology::END_WHERE_CLAUSE;
    if(searchType=="PDB")
    {
      query<<"VALUES ?pdb { \"" << searchTerm << "\" }\n";
    }
    else if(searchType=="Condensed_Sequence")
    {
      query<<"VALUES ?oligo_sequence { \"" << searchTerm << "\" }\n";
    }
    else if(searchType=="Oligo_Library")
    {
      // query<<"VALUES ?pdb { \"" << searchTerm << "\" }\n";
    }
    

    FormulateCURL(output_file_type, query.str());
}
// 
// PREFIX : <http://gmmo.uga.edu/#>
// PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
// PREFIX owl: <http://www.w3.org/2002/07/owl#>
// PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
// PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
// SELECT ?pdb ?title ?resolution ?Mean_B_Factor ?oligo_sequence ?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI
// WHERE
// {
// ?pdb_file :identifier ?pdb.
// ?pdb_file :hasTitle ?title.
// ?pdb_file :hasResolution ?resolution.
// ?pdb_file :hasAuthors ?authors.
// ?pdb_file :hasJournal ?journal.
// ?pdb_file :hasDOI ?DOI.
// ?pdb_file :hasPMID ?PMID.
// ?pdb_file :hasBFactor ?Mean_B_Factor.
// 
// ?oligo :oligoName ?oligo_sequence.
// ?oligo :oligoBFactor ?oligo_mean_B_Factor.
// }
// ORDER BY DESC(?resolution)