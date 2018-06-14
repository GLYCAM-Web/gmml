#include "../../../../includes/MolecularModeling/assembly.hpp"

std::string MolecularModeling::Assembly::QueryOntology(std::string searchType, std::string searchTerm, float resolution_min, float resolution_max, float b_factor_min, float b_factor_max, float oligo_b_factor_min, float oligo_b_factor_max, int page, int resultsPerPage, std::string sortBy, std::string url, std::string output_file_type)
{
    std::stringstream query;
    std::stringstream search;
    search << searchType;
       
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
    query << " REDUCED ?pdb ?oligo_sequence ?residue_links ?glycosidic_linkage ?title "
             "?resolution ?Mean_B_Factor ?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI "
             "(group_concat(distinct ?comment;separator=\"\\n\") as ?comments) "
             "(group_concat(distinct ?warning;separator=\"\\n\") as ?warnings) "
             "(group_concat(distinct ?error;separator=\"\\n\") as ?errors)\n";  
    query << Ontology::WHERE_CLAUSE;
    query << "?oligo        :oligoName              ?oligo_sequence.\n";
    if(search.str()=="Oligo_REGEX")
    {
      query << "FILTER regex(?oligo_sequence, \"" << searchTerm << "\")\n";    
    }
    if(search.str()=="Condensed_Sequence")
    {
      gmml::FindReplaceString(searchTerm, "[", "\\\\[");
      gmml::FindReplaceString(searchTerm, "]", "\\\\]");
      gmml::FindReplaceString(searchTerm, "-OH", "-ROH");
      query << "VALUES ?oligo_sequence { \"" << searchTerm << "\" }\n";
    }
    query << "?pdb_file     :identifier             ?pdb.\n";
    if(search.str()=="PDB")
    {
      query << "VALUES ?pdb { \"" << searchTerm << "\" }\n";
    }
    query << "?pdb_file     :hasTitle               ?title.\n";
    if(resolution_max == -1 && resolution_min == -1)
    {
      query << "OPTIONAL {";
    }
    query << "?pdb_file     :hasResolution          ?resolution.";
    if(resolution_max != -1)
    {
      query << "FILTER (" << resolution_max << " > ?resolution)";
    }
    if(resolution_min != -1)
    {
      query << "FILTER (" << resolution_min << " < ?resolution)";
    }
    if(resolution_max == -1 && resolution_min == -1)
    {
      query << "}";
    }
    query << "\n";
    query << "?pdb_file     :hasAuthors             ?authors.\n";
    query << "OPTIONAL {";
    query << "?pdb_file     :hasJournal             ?journal.}\n";
    query << "OPTIONAL {";
    query << "?pdb_file     :hasDOI                 ?DOI.}\n";
    query << "OPTIONAL {";
    query << "?pdb_file     :hasPMID                ?PMID.}\n";
    if(b_factor_max != -1)
    {
      query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.";
      query << "FILTER (" << b_factor_max << " > ?Mean_B_Factor)";
    }
    if(b_factor_min != -1)
    {
      query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.";
      query << "FILTER (" << b_factor_min << " < ?Mean_B_Factor)";
    }
    if(b_factor_max == -1 && b_factor_min == -1)
    {
      query << "OPTIONAL {";
      query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.}";
    }
    query << "\n";
    query << "OPTIONAL {";
    query << "?oligo        :oligoName              ?oligo_sequence.}\n";
    if(oligo_b_factor_max != -1)
    {
      query << "FILTER (" << oligo_b_factor_max << " > ?oligo_mean_B_Factor)";
    }
    if(oligo_b_factor_min != -1)
    {
      query << "FILTER (" << oligo_b_factor_min << " < ?oligo_mean_B_Factor)";
    }
    if(oligo_b_factor_max == -1 && oligo_b_factor_min == -1)
    {
      query << "OPTIONAL {";
      query << "?oligo        :oligoBFactor           ?oligo_mean_B_Factor.}";
    }
    query << "\n";
    query << "OPTIONAL {";
    query << "?pdb_file     :hasOligo	              ?oligo.}\n";
    query << "OPTIONAL {";
    query << "?oligo        :oligoResidueLinks      ?residue_links.}\n";
    query << "OPTIONAL {";
    query << "?linkage      :hasParent 	            ?oligo.}\n";
    query << "OPTIONAL {";
    query << "?linkage      :glycosidicLinkage      ?glycosidic_linkage.}\n";
    query << "OPTIONAL {";
    query << "?pdb_file      :hasNote    ?errorNote.\n";
    query << "?errornote	       :NoteType    \"error\".\n";
    query << "?errornote        :description ?error.}\n";
    query << "OPTIONAL {";
    query << "?pdb_file      :hasNote    ?warningNote.\n";
    query << "?warningNote	       :NoteType    \"warning\".\n";
    query << "?warningNote        :description ?warning.}\n";
    query << "OPTIONAL {";
    query << "?pdb_file      :hasNote    ?commentNote.\n";
    query << "?commentNote	       :NoteType    \"comment\".\n";
    query << "?commentNote        :description ?comment.}\n";
    query << Ontology::END_WHERE_CLAUSE << "\n";
    query << "ORDER BY  ?" << sortBy << "\n";
    if(resultsPerPage != -1)
    {
    query << "LIMIT  " << resultsPerPage << "\n";
    }
    query << "OFFSET " << resultsPerPage*(page - 1) << "\n";
    
    
  
    
    
    

    return FormulateCURLGF(output_file_type, query.str(), url);
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