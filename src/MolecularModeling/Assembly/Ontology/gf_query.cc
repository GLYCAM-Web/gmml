#include "../../../../includes/MolecularModeling/assembly.hpp"

std::string MolecularModeling::Assembly::QueryOntology(std::string searchType, std::string searchTerm, std::string url, std::string output_file_type)
{
    std::stringstream query;
    std::stringstream search;
    search << searchType;
    if(search.str()=="Oligo_REGEX")
    {      
      query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
      query << " REDUCED ?pdb ?oligo_sequence ?residue_links ?glycosidic_linkage ?title "
               "?resolution ?Mean_B_Factor ?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI "
               "(group_concat(distinct ?comment;separator=\"\\n\") as ?comments) "
               "(group_concat(distinct ?warning;separator=\"\\n\") as ?warnings) "
               "(group_concat(distinct ?error;separator=\"\\n\") as ?errors)\n";  
      query << Ontology::WHERE_CLAUSE;
      query << "?oligo        :oligoName              ?oligo_sequence.\n";
      
      query << "FILTER regex(?oligo_sequence, \"" << searchTerm << "\")\n";    
      query << "?pdb_file     :identifier             ?pdb.\n";
      query << "?pdb_file     :hasTitle               ?title.\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasResolution          ?resolution.}\n";
      query << "?pdb_file     :hasAuthors             ?authors.\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasJournal             ?journal.}\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasDOI                 ?DOI.}\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasPMID                ?PMID.}\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.}\n";
      query << "OPTIONAL {";
      query << "?oligo        :oligoName              ?oligo_sequence.}\n";
      query << "OPTIONAL {";
      query << "?oligo        :oligoBFactor           ?oligo_mean_B_Factor.}\n";
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
      query << Ontology::END_WHERE_CLAUSE;
    }
    else
    {
      query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
      query << " REDUCED ?pdb ?oligo_sequence ?residue_links ?glycosidic_linkage ?title "
               "?resolution ?Mean_B_Factor ?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI "
               "(group_concat(distinct ?comment;separator=\"\\n\") as ?comments) "
               "(group_concat(distinct ?warning;separator=\"\\n\") as ?warnings) "
               "(group_concat(distinct ?error;separator=\"\\n\") as ?errors)\n";  
      query << Ontology::WHERE_CLAUSE;
      query << "?oligo        :oligoName              ?oligo_sequence.\n";
      query << "?pdb_file     :identifier             ?pdb.\n";
      query << "?pdb_file     :hasTitle               ?title.\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasResolution          ?resolution.}\n";
      query << "?pdb_file     :hasAuthors             ?authors.\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasJournal             ?journal.}\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasDOI                 ?DOI.}\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasPMID                ?PMID.}\n";
      query << "OPTIONAL {";
      query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.}\n";
      query << "OPTIONAL {";
      query << "?oligo        :oligoName              ?oligo_sequence.}\n";
      query << "OPTIONAL {";
      query << "?oligo        :oligoBFactor           ?oligo_mean_B_Factor.}\n";
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
      query << Ontology::END_WHERE_CLAUSE;
      
      if(search.str()=="PDB")
      {
        query<<"VALUES ?pdb { \"" << searchTerm << "\" }\n";
      }
      else if(search.str()=="Condensed_Sequence")
      {
        gmml::FindReplaceString(searchTerm, "[", "\\\\[");
        gmml::FindReplaceString(searchTerm, "]", "\\\\]");
        gmml::FindReplaceString(searchTerm, "-OH", "-ROH");
        query<<"VALUES ?oligo_sequence { \"" << searchTerm << "\" }\n";
      }
    }
    
    
    

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