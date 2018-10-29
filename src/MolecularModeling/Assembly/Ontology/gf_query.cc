#include "../../../../includes/MolecularModeling/assembly.hpp"

std::string MolecularModeling::Assembly::MoreQuery(std::string pdb_id, std::string oligo_sequence, std::string oligo, std::string url, std::string output_file_type)
{
  std::stringstream query;
  
  query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
  query << " DISTINCT ?residue_links ?glycosidic_linkage ?title ?resolution ?Mean_B_Factor"
           "?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI"
           "(group_concat(distinct ?comment;separator=\"\\n\") as ?comments)" 
           "(group_concat(distinct ?warning;separator=\"\\n\") as ?warnings)" 
           "(group_concat(distinct ?error;separator=\"\\n\") as ?errors)";
  query << Ontology::WHERE_CLAUSE;
  query << "?pdb_file     :identifier    \"" << pdb_id << "\";\n";
  query << "              :hasOligo      ?oligo.\n";
  query << "FILTER regex(?oligo, \"" << oligo << "$\")\n";
  // gmml::FindReplaceString(oligo_sequence, "[", "\\\\[");
  // gmml::FindReplaceString(oligo_sequence, "]", "\\\\]");
  gmml::FindReplaceString(oligo_sequence, "-OH", "-ROH");
  query << "?oligo        :oligoName     \"" << oligo_sequence << "\".\n";
  query << "?pdb_file     :hasTitle               ?title;\n";
  query << "              :hasAuthors             ?authors.\n";
  query << "OPTIONAL {";
  query << "?pdb_file     :hasJournal             ?journal.}\n";
  query << "OPTIONAL {";
  query << "?pdb_file     :hasDOI                 ?DOI.}\n";
  query << "OPTIONAL {";
  query << "?pdb_file     :hasPMID                ?PMID.}\n";
  query << "OPTIONAL {";
  query << "?pdb_file     :hasResolution          ?resolution.}\n";
  query << "OPTIONAL {";
  query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.}\n";
  query << "OPTIONAL {";
  query << "?oligo        :oligoResidueLinks      ?residue_links.}\n";
  query << "OPTIONAL {";
  query << "?oligo        :oligoBFactor           ?oligo_mean_B_Factor.}\n";
  query << "OPTIONAL {";
  query << "?linkage      :hasParent 	            ?oligo;\n";
  query << "              :glycosidicLinkage      ?glycosidic_linkage.}\n";
  query << "OPTIONAL {";
  query << "?pdb_file       :hasNote       ?errorNote.\n";
  query << "?errorNote	    :NoteType      \"error\".\n";
  query << "?errorNote      :description   ?error.}\n";
  query << "OPTIONAL {";
  query << "?pdb_file       :hasNote       ?warningNote.\n";
  query << "?warningNote    :NoteType      \"warning\".\n";
  query << "?warningNote    :description   ?warning.}\n";
  query << "OPTIONAL {";
  query << "?pdb_file       :hasNote       ?commentNote.\n";
  query << "?commentNote    :NoteType      \"comment\".\n";
  query << "?commentNote    :description   ?comment.}\n";
  query << Ontology::END_WHERE_CLAUSE << "\n";
  
  
  return FormulateCURLGF(output_file_type, query.str(), url);

}


std::string MolecularModeling::Assembly::QueryOntology(std::string searchType, std::string searchTerm, float resolution_min, float resolution_max, float b_factor_min, float b_factor_max, float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning, int isComment, int page, int resultsPerPage, std::string sortBy, std::string url, std::string output_file_type)
{
    std::stringstream query;
    std::stringstream search;
    search << searchType;
       
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
    query << " DISTINCT ?pdb ?oligo ?oligo_sequence \n";
    if(isComment == 1)
    {
      query << "(group_concat(distinct ?comment;separator=\"\\n\") as ?comments) ";
    }
    if(isWarning == 1)
    {
      query << "(group_concat(distinct ?warning;separator=\"\\n\") as ?warnings) ";
    }
    if(isError == 1)
    {
       query << "(group_concat(distinct ?error;separator=\"\\n\") as ?errors)\n";  
    }     
    query << Ontology::WHERE_CLAUSE;
    query << "?pdb_file     :identifier             ?pdb.\n";
    if(search.str()=="PDB")
    {
      query << "VALUES ?pdb { \"" << searchTerm << "\" }\n";
    }
    if((resolution_max != -1) | (resolution_min != -1))
    {
      query << "?pdb_file     :hasResolution          ?resolution.\n";
    }
    if(resolution_max != -1)
    {
      query << "FILTER (" << resolution_max << " > ?resolution)\n";
    }
    if(resolution_min != -1)
    {
      query << "FILTER (" << resolution_min << " < ?resolution)\n";
    }
    if((b_factor_max != -1) | (b_factor_min != -1))
    {
      query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.\n";
    }
    if(b_factor_max != -1)
    {
      query << "FILTER (" << b_factor_max << " > ?Mean_B_Factor)\n";
    }
    if(b_factor_min != -1)
    {
      query << "FILTER (" << b_factor_min << " < ?Mean_B_Factor)\n";
    }
    query << "?pdb_file     :hasOligo               ?oligo.\n";
    query << "?oligo        :oligoName              ?oligo_sequence.\n";
    if(search.str()=="Oligo_REGEX")
    {
      gmml::FindReplaceString(searchTerm, "[", "\\\\[");
      gmml::FindReplaceString(searchTerm, "]", "\\\\]");
      gmml::FindReplaceString(searchTerm, "-OH", "-ROH");
      query << "FILTER regex(?oligo_sequence, \"" << searchTerm << "\")\n";    
    }
    if(search.str()=="Condensed_Sequence")
    {
      gmml::FindReplaceString(searchTerm, "[", "\\\\[");
      gmml::FindReplaceString(searchTerm, "]", "\\\\]");
      gmml::FindReplaceString(searchTerm, "-OH", "-ROH");
      query << "VALUES ?oligo_sequence { \"" << searchTerm << "\" }\n";
    }
    if((oligo_b_factor_max != -1) | (oligo_b_factor_min != -1))
    {
      query << "?oligo        :oligoBFactor           ?oligo_mean_B_Factor.\n";
    }
    if(oligo_b_factor_max != -1)
    {
      query << "FILTER (" << oligo_b_factor_max << " > ?oligo_mean_B_Factor)\n";
    }
    if(oligo_b_factor_min != -1)
    {
      query << "FILTER (" << oligo_b_factor_min << " < ?oligo_mean_B_Factor)\n";
    }
    if(isError == 1)
    {
      query << "?pdb_file       :hasNote       ?errorNote.\n";
      query << "?errorNote	    :NoteType      \"error\".\n";
      query << "?errorNote      :description   ?error.\n";
    }
    if(isWarning == 1)
    {
      query << "?pdb_file       :hasNote       ?warningNote.\n";
      query << "?warningNote    :NoteType      \"warning\".\n";
      query << "?warningNote    :description   ?warning.\n";
    }
    if(isComment == 1)
    {
      query << "?pdb_file       :hasNote       ?commentNote.\n";
      query << "?commentNote    :NoteType      \"comment\".\n";
      query << "?commentNote    :description   ?comment.\n";
    }
  
    query << Ontology::END_WHERE_CLAUSE << "\n";
    query << "ORDER BY  ?" << sortBy << "\n";
    if(resultsPerPage != -1)
    {
    query << "LIMIT  " << resultsPerPage << "\n";
    }
    query << "OFFSET " << resultsPerPage*(page - 1) << "\n";
    
    
    return FormulateCURLGF(output_file_type, query.str(), url);
}

std::string MolecularModeling::Assembly::ontologyPDBDownload(std::string searchType, std::string searchTerm, float resolution_min, float resolution_max, float b_factor_min, float b_factor_max, float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning, int isComment, std::string sortBy, std::string url, std::string output_file_type)
{
  std::stringstream query;
  std::stringstream search;
  search << searchType;
     
  query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
  query << " DISTINCT ?PDB_ID \n";
  query << "(group_concat(distinct ?oligo_sequence;separator=\"\\n\") as ?Oligosaccharides) ";
  if(isComment == 1)
  {
    query << "(group_concat(distinct ?comment;separator=\"\\n\") as ?comments) ";
  }
  if(isWarning == 1)
  {
    query << "(group_concat(distinct ?warning;separator=\"\\n\") as ?warnings) ";
  }
  if(isError == 1)
  {
     query << "(group_concat(distinct ?error;separator=\"\\n\") as ?errors)\n";  
  }     
  query << Ontology::WHERE_CLAUSE;
  query << "?pdb_file     :identifier             ?PDB_ID.\n";
  if(search.str()=="PDB")
  {
    query << "VALUES ?PDB_ID { \"" << searchTerm << "\" }\n";
  }
  if((resolution_max != -1) | (resolution_min != -1))
  {
    query << "?pdb_file     :hasResolution          ?resolution.\n";
  }
  if(resolution_max != -1)
  {
    query << "FILTER (" << resolution_max << " > ?resolution)\n";
  }
  if(resolution_min != -1)
  {
    query << "FILTER (" << resolution_min << " < ?resolution)\n";
  }
  if((b_factor_max != -1) | (b_factor_min != -1))
  {
    query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.\n";
  }
  if(b_factor_max != -1)
  {
    query << "FILTER (" << b_factor_max << " > ?Mean_B_Factor)\n";
  }
  if(b_factor_min != -1)
  {
    query << "FILTER (" << b_factor_min << " < ?Mean_B_Factor)\n";
  }
  query << "?pdb_file     :hasOligo               ?oligo.\n";
  query << "?oligo        :oligoName              ?oligo_sequence.\n";
  if(search.str()=="Oligo_REGEX")
  {
    gmml::FindReplaceString(searchTerm, "[", "\\\\[");
    gmml::FindReplaceString(searchTerm, "]", "\\\\]");
    gmml::FindReplaceString(searchTerm, "-OH", "-ROH");
    query << "FILTER regex(?oligo_sequence, \"" << searchTerm << "\")\n";    
  }
  if(search.str()=="Condensed_Sequence")
  {
    gmml::FindReplaceString(searchTerm, "[", "\\\\[");
    gmml::FindReplaceString(searchTerm, "]", "\\\\]");
    gmml::FindReplaceString(searchTerm, "-OH", "-ROH");
    query << "VALUES ?oligo_sequence { \"" << searchTerm << "\" }\n";
  }
  if((oligo_b_factor_max != -1) | (oligo_b_factor_min != -1))
  {
    query << "?oligo        :oligoBFactor           ?oligo_mean_B_Factor.\n";
  }
  if(oligo_b_factor_max != -1)
  {
    query << "FILTER (" << oligo_b_factor_max << " > ?oligo_mean_B_Factor)\n";
  }
  if(oligo_b_factor_min != -1)
  {
    query << "FILTER (" << oligo_b_factor_min << " < ?oligo_mean_B_Factor)\n";
  }
  if(isError == 1)
  {
    query << "?pdb_file       :hasNote       ?errorNote.\n";
    query << "?errorNote	    :NoteType      \"error\".\n";
    query << "?errorNote      :description   ?error.\n";
  }
  if(isWarning == 1)
  {
    query << "?pdb_file       :hasNote       ?warningNote.\n";
    query << "?warningNote    :NoteType      \"warning\".\n";
    query << "?warningNote    :description   ?warning.\n";
  }
  if(isComment == 1)
  {
    query << "?pdb_file       :hasNote       ?commentNote.\n";
    query << "?commentNote    :NoteType      \"comment\".\n";
    query << "?commentNote    :description   ?comment.\n";
  }

  query << Ontology::END_WHERE_CLAUSE << "\n";
  query << "ORDER BY  ?" << sortBy << "\n";
    
  return FormulateCURLGF(output_file_type, query.str(), url);
}

std::string MolecularModeling::Assembly::ontologyDownload(std::string searchType, std::string searchTerm, float resolution_min, float resolution_max, float b_factor_min, float b_factor_max, float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning, int isComment, std::string sortBy, std::string url, std::string output_file_type)
{
  std::stringstream query;
  std::stringstream search;
  search << searchType;

  query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
  query << " DISTINCT ?pdb ?oligo_sequence ?residue_links ?title ?resolution ?Mean_B_Factor"
           "?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI"
           "(group_concat(distinct ?comment;separator=\"\\n\") as ?comments)"
           "(group_concat(distinct ?warning;separator=\"\\n\") as ?warnings)"
           "(group_concat(distinct ?error;separator=\"\\n\") as ?errors)";
           query << Ontology::WHERE_CLAUSE;
           query << "?pdb_file     :identifier             ?pdb.\n";
           query << "?pdb_file     :hasTitle               ?title;\n";
           query << "              :hasAuthors             ?authors.\n";
           query << "OPTIONAL {";
           query << "?pdb_file     :hasJournal             ?journal.}\n";
           query << "OPTIONAL {";
           query << "?pdb_file     :hasDOI                 ?DOI.}\n";
           query << "OPTIONAL {";
           query << "?pdb_file     :hasPMID                ?PMID.}\n";
           if(search.str()=="PDB")
           {
             query << "VALUES ?pdb { \"" << searchTerm << "\" }\n";
           }
           query << "OPTIONAL {";
           query << "?pdb_file     :hasResolution          ?resolution.\n";
           if(resolution_max != -1)
           {
             query << "FILTER (" << resolution_max << " > ?resolution)\n";
           }
           if(resolution_min != -1)
           {
             query << "FILTER (" << resolution_min << " < ?resolution)\n";
           }
           query << "}\n";
           query << "OPTIONAL {";
           query << "?pdb_file     :hasBFactor             ?Mean_B_Factor.\n";
           if(b_factor_max != -1)
           {
             query << "FILTER (" << b_factor_max << " > ?Mean_B_Factor)\n";
           }
           if(b_factor_min != -1)
           {
             query << "FILTER (" << b_factor_min << " < ?Mean_B_Factor)\n";
           }
           query << "}\n";
           query << "?pdb_file     :hasOligo               ?oligo.\n";
           query << "?oligo        :oligoName              ?oligo_sequence.\n";
           if(search.str()=="Oligo_REGEX")
           {
             gmml::FindReplaceString(searchTerm, "[", "\\\\[");
             gmml::FindReplaceString(searchTerm, "]", "\\\\]");
             gmml::FindReplaceString(searchTerm, "-OH", "-ROH");
             query << "FILTER regex(?oligo_sequence, \"" << searchTerm << "\")\n";
           }
           if(search.str()=="Condensed_Sequence")
           {
             gmml::FindReplaceString(searchTerm, "[", "\\\\[");
             gmml::FindReplaceString(searchTerm, "]", "\\\\]");
             gmml::FindReplaceString(searchTerm, "-OH", "-ROH");
             query << "VALUES ?oligo_sequence { \"" << searchTerm << "\" }\n";
           }
           query << "OPTIONAL {";
           query << "?oligo        :oligoResidueLinks      ?residue_links.}\n";
           query << "OPTIONAL {";
           query << "?oligo        :oligoBFactor           ?oligo_mean_B_Factor.\n";
           if(oligo_b_factor_max != -1)
           {
             query << "FILTER (" << oligo_b_factor_max << " > ?oligo_mean_B_Factor)\n";
           }
           if(oligo_b_factor_min != -1)
           {
             query << "FILTER (" << oligo_b_factor_min << " < ?oligo_mean_B_Factor)\n";
           }
           query << "}\n";
           if(isError != 1)
           {
             query << "OPTIONAL {";
           }
           query << "?pdb_file       :hasNote       ?errorNote.\n";
           query << "?errorNote	    :NoteType      \"error\".\n";
           query << "?errorNote      :description   ?error.\n";
           if(isError != 1)
           {
             query << "}\n";
           }
           if(isWarning != 1)
           {
             query << "OPTIONAL {";
           }
           query << "?pdb_file       :hasNote       ?warningNote.\n";
           query << "?warningNote    :NoteType      \"warning\".\n";
           query << "?warningNote    :description   ?warning.\n";
           if(isWarning != 1)
           {
             query << "}\n";
           }
           if(isComment != 1)
           {
             query << "OPTIONAL {";
           }
           query << "?pdb_file       :hasNote       ?commentNote.\n";
           query << "?commentNote    :NoteType      \"comment\".\n";
           query << "?commentNote    :description   ?comment.\n";
           if(isComment != 1)
           {
             query << "}\n";
           }

           query << Ontology::END_WHERE_CLAUSE << "\n";
           query << "ORDER BY  ?" << sortBy << "\n";
           





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