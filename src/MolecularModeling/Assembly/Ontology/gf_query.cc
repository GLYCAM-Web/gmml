#include "../../../../includes/MolecularModeling/assembly.hpp"

//For an example query with some explaination, see the bottom of this file.  For sparql query information, see https://www.w3.org/TR/rdf-sparql-query/ (It is not the greatest documentation but it helps)


std::string MolecularModeling::Assembly::MoreQuery(std::string pdb_id, std::string oligo_sequence, std::string oligo, std::string url, std::string output_file_type)
{ // This function runs a full query on a single result, which is unique given the pdb_id, oligo_sequence, and oligo (which is the oligosaccharide number in the PDB in case it has identical sugars found)

  std::stringstream query;
  query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
  query << " DISTINCT ?residue_links" /*?glycosidic_linkage*/ "?title ?resolution ?Mean_B_Factor"
           "?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI ?pdb_coordinates ?ProteinID"
           "(group_concat(distinct ?comment;separator=\"\\n\") as ?comments)"
           "(group_concat(distinct ?warning;separator=\"\\n\") as ?warnings)"
           "(group_concat(distinct ?error;separator=\"\\n\") as ?errors)";
  query << Ontology::WHERE_CLAUSE;
  query << "?pdb_file     :identifier    \"" << pdb_id << "\";\n";
  query << "              :hasOligo      ?oligo.\n";
  query << "FILTER regex(?oligo, \"" << oligo << "$\")\n";
  gmml::FindReplaceString(oligo_sequence, "-OH", "-ROH");
  query << "?oligo        :oligoIUPACname     \"" << oligo_sequence << "\".\n";
  query << "?pdb_file     :hasTitle               ?title;\n";
  query << "              :hasAuthors             ?authors.\n";
  query << "OPTIONAL {";
  query << "?pdb_file     :hasJournal             ?journal.}\n";
  query << "OPTIONAL {";
  query << "?pdb_file     :hasProteinID           ?ProteinID.}\n";
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
  query << "?oligo        :PDBfile           ?pdb_coordinates.\n";
  query << "?oligo        :hasMono            ?mono.\n";
  // query << "OPTIONAL {";
  // query << "?linkage      :hasParent 	            ?oligo;\n";
  // query << "              :glycosidicLinkage      ?glycosidic_linkage.}\n";
  query << "OPTIONAL {";
  query << "?mono       :hasNote       ?errorNote.\n";
  query << "?errorNote	    :NoteType      \"error\".\n";
  query << "?errorNote      :description   ?error.}\n";
  query << "OPTIONAL {";
  query << "?mono       :hasNote       ?warningNote.\n";
  query << "?warningNote    :NoteType      \"warning\".\n";
  query << "?warningNote    :description   ?warning.}\n";
  query << "OPTIONAL {";
  query << "?mono       :hasNote       ?commentNote.\n";
  query << "?commentNote    :NoteType      \"comment\".\n";
  query << "?commentNote    :description   ?comment.}\n";
  //add info for coordinates here
  query << Ontology::END_WHERE_CLAUSE << "\n";

  // gmml::log(__LINE__, __FILE__, gmml::INF, query.str());
  // std::cout << "\n" << query.str() << "\n";
  return FormulateCURLGF(output_file_type, query.str(), url);

}

// std::string MolecularModeling::Assembly::BranchQuery(std::string searchType, std::string searchTerm, float resolution_min, float resolution_max, float b_factor_min, float b_factor_max, float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning, int isComment, int page, int resultsPerPage, std::string sortBy, std::string url, std::string output_file_type)
// {
//   //TODO figure out how to split up oligo sequence between monos handling branching brackets
// }

std::string MolecularModeling::Assembly::QueryOntology(std::string searchType, std::string searchTerm, float resolution_min, float resolution_max, float b_factor_min, float b_factor_max, float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning, int isComment, int isLigand, int isGlycomimetic, int isNucleotide, std::string aglycon, std::string count, int page, int resultsPerPage, std::string sortBy, std::string url, std::string output_file_type)
{   //This function runs a basic query, looking only for ?pdb (PDB_ID), ?oligo (Oligosaccharides are assigned numbers when they are found, ie oligo_1),
    //and ?oligo_sequence (Condensed sequence).  These three variables together are unique for each result.  This function also takes in all of the possible
    //filter variables to return filtered results when updating via ajax

    std::stringstream query;
    std::stringstream search;
    search << searchType;

    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE;
    if(count == "TRUE")
    {
      query << " DISTINCT (COUNT(?oligo) as ?count) \n";
    }
    else
    {
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
    query << "?oligo        :hasMono                ?mono.\n";
    query << "?oligo        :oligoIUPACname              ?oligo_sequence.\n";
    query << "FILTER (!regex(?oligo_sequence, \"" << "Unknown" << "\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"" << "HEM" << "\"))\n";
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

    if(isLigand == 1)
    {
      query << "FILTER (!regex(?oligo_sequence, \"-ASN$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-THR$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-SER$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-LYZ$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-HYP$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-TYR$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-CYS$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-TRP$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-LYS$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-HIS$\"))\n";
    }
    else if(isLigand == 0)
    {
      query << "FILTER (!regex(?oligo_sequence, \"-ROH$\"))\n";
      query << "FILTER (!regex(?oligo_sequence, \"-OME$\"))\n";
      query << "?oligo    :oligoSequenceName     ?sequenceName.\n";
      query << "FILTER (!regex(?sequenceName, \"-Unknown$\"))\n";
    }
    if(isNucleotide == 1)
    {
      query << "?mono         :isNucleotide  \"true\"\n";
    }
    else if(isNucleotide == 0)
    {
      query << "?mono         :isNucleotide  \"false\"\n";
    }
    if(isGlycomimetic == 1)
    {
      query << "FILTER regex(?oligo_sequence, \"<R\")\n";
    }
    else if(isGlycomimetic == 0)
    {
      query << "FILTER (!regex(?oligo_sequence, \"<R\"))\n";
    }

    if(aglycon.length() > 0)
    {
      query << "FILTER regex(?oligo_sequence, \"" << aglycon << "$\")\n";
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
      query << "?mono       :hasNote       ?errorNote.\n";
      query << "?errorNote	    :NoteType      \"error\".\n";
      query << "?errorNote      :description   ?error.\n";
    }
    else if(isError == 0)
    {
      query << "OPTIONAL {";
      query << "?mono       :hasNote       ?errorNote.}\n";
      query << "FILTER NOT EXISTS { ?errorNote :NoteType \"error\".}\n";
    }
    if(isWarning == 1)
    {
      query << "?mono       :hasNote       ?warningNote.\n";
      query << "?warningNote    :NoteType      \"warning\".\n";
      query << "?warningNote    :description   ?warning.\n";
    }
    else if(isWarning == 0)
    {
      query << "OPTIONAL {";
      query << "?mono       :hasNote       ?warningNote.}\n";
      query << "FILTER NOT EXISTS { ?warningNote :NoteType \"warning\".}\n";
    }
    if(isComment == 1)
    {
      query << "?mono       :hasNote       ?commentNote.\n";
      query << "?commentNote    :NoteType      \"comment\".\n";
      query << "?commentNote    :description   ?comment.\n";
    }
    else if(isComment == 0)
    {
      query << "OPTIONAL {";
      query << "?mono       :hasNote       ?commentNote.}\n";
      query << "FILTER NOT EXISTS { ?commentNote :NoteType \"comment\".}\n";
    }

    query << Ontology::END_WHERE_CLAUSE << "\n";
    if(count != "TRUE")
    {
      query << "ORDER BY  ?" << sortBy << "\n";
      if(resultsPerPage != -1)
      {
      query << "LIMIT  " << resultsPerPage << "\n";
      }
      query << "OFFSET " << resultsPerPage*(page - 1) << "\n";
    }
    // gmml::log(__LINE__, __FILE__, gmml::INF, query.str());

    return FormulateCURLGF(output_file_type, query.str(), url);
}

std::string MolecularModeling::Assembly::ontologyPDBDownload(std::string searchType, std::string searchTerm, float resolution_min, float resolution_max, float b_factor_min, float b_factor_max, float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning, int isComment, int isLigand, int isGlycomimetic, int isNucleotide, std::string aglycon, std::string count, int page, int resultsPerPage, std::string sortBy, std::string url, std::string output_file_type)
{ // This query creates a list of unique PDB_IDs given all of the user specified filters, and returns a CSV which is downloaded
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
  if(isLigand == 1)
  {
    query << "FILTER (!regex(?oligo_sequence, \"-ASN$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-THR$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-SER$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-LYZ$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-HYP$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-TYR$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-CYS$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-TRP$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-LYS$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-HIS$\"))\n";
  }
  else if(isLigand == 0)
  {
    query << "FILTER (!regex(?oligo_sequence, \"-ROH$\"))\n";
    query << "FILTER (!regex(?oligo_sequence, \"-OME$\"))\n";
    query << "?oligo    :oligoSequenceName     ?sequenceName.\n";
    query << "FILTER (!regex(?sequenceName, \"-Unknown$\"))\n";
  }
  if(isNucleotide == 1)
  {
    query << "?mono         :isNucleotide  \"true\"\n";
  }
  else if(isNucleotide == 0)
  {
    query << "?mono         :isNucleotide  \"false\"\n";
  }
  if(isGlycomimetic == 1)
  {
    query << "FILTER regex(?oligo_sequence, \"<R\")\n";
  }
  else if(isGlycomimetic == 0)
  {
    query << "FILTER (!regex(?oligo_sequence, \"<R\"))\n";
  }

  if(aglycon.length() > 0)
  {
    query << "FILTER regex(?oligo_sequence, \"" << aglycon << "$\")\n";
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

std::string MolecularModeling::Assembly::ontologyDownload(std::string searchType, std::string searchTerm, float resolution_min, float resolution_max, float b_factor_min, float b_factor_max, float oligo_b_factor_min, float oligo_b_factor_max, int isError, int isWarning, int isComment, int isLigand, int isGlycomimetic, int isNucleotide, std::string aglycon, std::string count, int page, int resultsPerPage, std::string sortBy, std::string url, std::string output_file_type)
{ //This is a complete (and therefore slow) query that is a combination of moreQuery() and QueryOntology().  It filters the database by user input, and returns a CSV with all of the data for download.
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
           if(isLigand == 1)
           {
             query << "FILTER (!regex(?oligo_sequence, \"-ASN$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-THR$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-SER$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-LYZ$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-HYP$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-TYR$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-CYS$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-TRP$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-LYS$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-HIS$\"))\n";
           }
           else if(isLigand == 0)
           {
             query << "FILTER (!regex(?oligo_sequence, \"-ROH$\"))\n";
             query << "FILTER (!regex(?oligo_sequence, \"-OME$\"))\n";
             query << "?oligo    :oligoSequenceName     ?sequenceName.\n";
             query << "FILTER (!regex(?sequenceName, \"-Unknown$\"))\n";
           }
           if(isNucleotide == 1)
           {
             query << "?mono         :isNucleotide  \"true\"\n";
           }
           else if(isNucleotide == 0)
           {
             query << "?mono         :isNucleotide  \"false\"\n";
           }
           if(isGlycomimetic == 1)
           {
             query << "FILTER regex(?oligo_sequence, \"<R\")\n";
           }
           else if(isGlycomimetic == 0)
           {
             query << "FILTER (!regex(?oligo_sequence, \"<R\"))\n";
           }

           if(aglycon.length() > 0)
           {
             query << "FILTER regex(?oligo_sequence, \"" << aglycon << "$\")\n";
           }
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

// To test new queries, you should go to your dev site's online virtuoso query.
// To find the IP address of your Virtuoso database, run docker inspect (YOUR_USERNAME)_gw_virt | grep IPAddress
// Then go to that IP address, followed by :8890/sparql.  SO for me that's http://172.16.3.8:8890/sparql
// To query our ontology, you need to use the prefixes below.  This tells the ontology what the vocabulary means.
// SELECT is all of the data you want to pull out, and each variable must be present in the WHERE{} part of the query.
// Below is a sample query, and above if you follow the code it generates a couple other queries.
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
// FILTER (!regex(?oligo_IUPAC, "- Unknown$"))
// ?oligo :oligoName ?oligo_sequence.
// ?oligo :oligoBFactor ?oligo_mean_B_Factor.
// }
// ORDER BY DESC(?resolution)


//Here's another that I ran for Rob to get all PDBs with non furanose (!regex line) sugars with unercognized side chains (symbolized as <R)
//
// PREFIX : <http://gmmo.uga.edu/#>
// PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
// PREFIX owl: <http://www.w3.org/2002/07/owl#>
// PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
// PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
// SELECT DISTINCT ?pdb ?oligo (group_concat(distinct ?oligo_sequence;separator="\n") as ?Oligosaccharides)
// WHERE
// {
// ?pdb_file :identifier ?pdb.
// ?pdb_file :hasOligo ?oligo.
// ?oligo :oligoName ?oligo_sequence.
// FILTER regex(?oligo_sequence, ".*<R.*")
// FILTER (!regex(?oligo_sequence, ".*f.*"))
// }

//More filters that I am saving here just in case
// ?oligo :oligoIUPACname ?oligo_IUPAC.
// FILTER regex(?oligo_sequence, ".*-Unknown$")
// FILTER (!regex(?oligo_sequence, ".*<R.*"))
// FILTER (!regex(?oligo_IUPAC, "- Unknown$"))
// FILTER (!regex(?oligo_IUPAC, "-$"))
// FILTER (!regex(?oligo_IUPAC, "-ASN$"))
// FILTER (!regex(?oligo_IUPAC, "-THR$"))
// FILTER (!regex(?oligo_IUPAC, "-SER$"))


// PREFIX : <http://gmmo.uga.edu/#>
// PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
// PREFIX owl: <http://www.w3.org/2002/07/owl#>
// PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
// PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
// SELECT ?pdb ?oligo_sequence
// WHERE
// {
// ?pdb_file :identifier ?pdb.
// ?pdb_file :hasOligo ?oligo.
// ?oligo :oligoName ?oligo_sequence.

// SELECT DISTINCT ?residue_links?title ?resolution ?Mean_B_Factor?oligo_mean_B_Factor ?authors ?journal ?PMID ?DOI ?pdb_coordinates ?ProteinID(group_concat(distinct ?comment;separator="\n") as ?comments)(group_concat(distinct ?warning;separator="\n") as ?warnings)(group_concat(distinct ?error;separator="\n") as ?errors)WHERE {
// ?pdb_file     :identifier    "1A14";
//               :hasOligo      ?oligo.
// FILTER regex(?oligo, "oligo1$")
// ?oligo        :oligoIUPACname     "DManpa1-2DManpa1-2DManpa1-3DManpb1-4DGlcpNAcb1-4DGlcpNAcb1-ASN".
// ?pdb_file     :hasTitle               ?title;
//               :hasAuthors             ?authors.
// OPTIONAL {?pdb_file     :hasJournal             ?journal.}
// OPTIONAL {?pdb_file     :hasProteinID           ?ProteinID.}
// OPTIONAL {?pdb_file     :hasDOI                 ?DOI.}
// OPTIONAL {?pdb_file     :hasPMID                ?PMID.}
// OPTIONAL {?pdb_file     :hasResolution          ?resolution.}
// OPTIONAL {?pdb_file     :hasBFactor             ?Mean_B_Factor.}
// OPTIONAL {?oligo        :oligoResidueLinks      ?residue_links.}
// OPTIONAL {?oligo        :oligoBFactor           ?oligo_mean_B_Factor.}
// ?oligo        :PDBfile           ?pdb_coordinates.
// ?oligo        :hasMono            ?mono.
// OPTIONAL {?mono       :hasNote       ?errorNote.
// ?errorNote	    :NoteType      "error".
// ?errorNote      :description   ?error.}
// OPTIONAL {?mono       :hasNote       ?warningNote.
// ?warningNote    :NoteType      "warning".
// ?warningNote    :description   ?warning.}
// OPTIONAL {?mono       :hasNote       ?commentNote.
// ?commentNote    :NoteType      "comment".
// ?commentNote    :description   ?comment.}
// }

//New query format for branched oligo searching
//Looking for DGlcpNAcb1-4[LFucpa1-3]DGlcpNAcb
// PREFIX : <http://gmmo.uga.edu/#>
// PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
// PREFIX owl: <http://www.w3.org/2002/07/owl#>
// PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
// PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
// SELECT DISTINCT ?pdb (group_concat(distinct ?Oligosaccharide;separator="\n") as ?Oligosaccharides)
// WHERE {
// ?pdb_file     :identifier    ?pdb;
//               :hasOligo      ?oligo.
// ?pdb_file     :hasTitle               ?title;
//               :hasAuthors             ?authors.
// ?oligo        :oligoIUPACname     ?Oligosaccharide.
// ?oligo        :hasMono            ?mono;
//               :hasSequenceResidue ?residue1.
// ?residue1     :monosaccharideShortName  """LFucpa""".
// 
// ?residue1     :is1-3ConnectedTo      ?residue2.
// ?residue2     :monosaccharideShortName """DGlcpNAcb""".
// ?residue3 :is1-4ConnectedTo ?residue2.
// ?residue3 :monosaccharideShortName """DGlcpNAcb""".
// }
