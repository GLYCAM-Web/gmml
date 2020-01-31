#include "../../../../includes/MolecularModeling/assembly.hpp"
#include "../../../../includes/MolecularModeling/Graph/graph.hpp"

#include <regex>

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

GraphDS::Graph MolecularModeling::Assembly::CreateQueryStringGraph(std::string queryString)
{
  int local_debug = 1;
  std::cout << "Starting graph creation\n";
  GraphDS::Graph graph;
  std::vector<parsedString> parsedVector;
  std::vector<std::string> nodeStrings;
  //split string into mono, linkage, and brackets
  //IE DGlcpNAcb1-4[LFucpa1-4]DGlcpNAcb1-ASN is split into:
  //{DGlcpNAcb, 1-4,     [,    LFucpa, 1-4,     ],    DGlcpNAcb, 1-,      ASN} (labels)
  //{Node1,     NULL,    NULL, Node2,  NULL,    NULL, Node3,     NULL,    Node4} (nodes)
  //{NULL,      Edge1-3, NULL, NULL,   Edge2-3, NULL, NULL,      Edge3-4, NULL} (edges)
  std::string labelStr;
  std::size_t linkageStart, linkageEnd, dashLocation;
  dashLocation = queryString.find("-");
  if(isdigit(queryString[dashLocation - 1]))
  {
    linkageStart = dashLocation - 1;
  }
  else
  {
    linkageStart = dashLocation;
  }
  if(isdigit(queryString[dashLocation + 1]))
  {
    linkageEnd = dashLocation + 1;
  }
  else
  {
    linkageEnd = dashLocation;
  }
  
  //Fill Vector, add nodes & labels
  for(unsigned int i=0; i<queryString.size(); i++)
  {
    // std::cout << i << " : " << queryString[i] << "\n";
    if((queryString[i] == '[')||(queryString[i] == ']'))
    {
      // std::cout << "There's a bracket at: " << i << "\n";
      labelStr.push_back(queryString[i]);
      parsedVector.push_back(parsedString(labelStr));
      labelStr = "";
      continue;
    }
    else if(i == linkageStart)
    {
      while(i +1 != linkageEnd + 1)
      {
        // std::cout << i << " : " << linkageEnd+1 << "\n";
        labelStr.push_back(queryString[i]);
        i++;
      }
      labelStr.push_back(queryString[i]);
      GraphDS::Edge* thisEdge = new GraphDS::Edge;
      thisEdge->AddEdgeLabel(labelStr);
      parsedVector.push_back(parsedString(labelStr, thisEdge));
      labelStr = "";
      dashLocation = queryString.find("-", linkageEnd + 1);
      if(dashLocation != std::string::npos)
      {
        if(isdigit(queryString[dashLocation - 1]))
        {
          linkageStart = dashLocation - 1;
        }
        else
        {
          linkageStart = dashLocation;
        }
        if(isdigit(queryString[dashLocation + 1]))
        {
          linkageEnd = dashLocation + 1;
        }
        else
        {
          linkageEnd = dashLocation;
        }
      }
      else
      {
        dashLocation = -1;
        linkageStart = dashLocation;
        linkageEnd = dashLocation;
      }
    }
    else
    {
      if(linkageStart != -1)
      {
        while(i+1 != linkageStart)
        {
          labelStr.push_back(queryString[i]);
          i++;
        }
        labelStr.push_back(queryString[i]);
        GraphDS::Node* newNode = new GraphDS::Node;
        nodeStrings.push_back(labelStr);
        void* ptr = &(nodeStrings[nodeStrings.size() -1][0]);
        newNode->SetNodeValue(ptr);
        newNode->SetNodeId(labelStr);
        parsedVector.push_back(parsedString(labelStr, newNode));
        labelStr = "";
      }
      else
      {
        while(i < queryString.size())
        {
          labelStr.push_back(queryString[i]);
          i++;
        }
        GraphDS::Node* newNode = new GraphDS::Node;
        nodeStrings.push_back(labelStr);
        void* ptr = &(nodeStrings[nodeStrings.size() -1][0]);
        newNode->SetNodeValue(ptr);
        newNode->SetNodeId(labelStr);
        parsedVector.push_back(parsedString(labelStr, newNode));
        labelStr = "";
      }
    }
  }
  
  //Go through vector, add edges b/w nodes & add nodes & edges to graph
  ConnectNodes(0, static_cast<int>(parsedVector.size()), parsedVector, graph);
  
  
  
  //Trying something different
  // 
  // //Creating a char* vector so I can remove the branches as they are made into nodes without destroying the 
  // //original string.  It also allows pointers to residues that haven't been turned into nodes yet
  // //Vector.erase(start, end) will remove the bits from start to end, in this case the branch
  // 
  // std::vector<char*> qStringPointer;
  // for(std::string::iterator it = queryString.begin(); it != queryString.end(); it++)
  // {
  //   char* c = &(*it);
  //   qStringPointer.push_back(c);
  // }
  // if(local_debug > 0)
  // {
  //   for(unsigned int i = 0; i < qStringPointer.size(); i++)
  //   {
  //     std::cout << *qStringPointer[i];
  //   }
  //   std::cout << "\n";
  // }
  // std::regex branchEndPattern("([1-2]-[1-9])(])");//will match a bracket after a linkage, so the end bracket of a branch
  // std::smatch matchLocation;
  // std::size_t numStartBrackets = std::count(queryString.begin(), queryString.end(), '[');
  // std::size_t numEndBrackets = std::count(queryString.begin(), queryString.end(), ']');
  // if(numStartBrackets != numEndBrackets)
  // {
  //   std::cerr << "This function will get stuck in an infinite loop of there are more or less open brackets than close brackets.\n";
  //   std::cerr << "Exiting CreateQueryStringGraph() now\n";
  //   return graph;
  // }
  // std::string queryStringCopy = queryString;
  // while((queryStringCopy.size() != 0) && (numStartBrackets == numEndBrackets))
  // {  
  //   // Below will return the first match in matchLocation.str(0)
  //   // Each pattern in () in the regex object will be returned in the following elements of the array
  //   // In this case, the linkage will be in matchLocation.str(1) and the bracket will be matchLocation.str(2)
  //   std::regex_search(queryStringCopy, matchLocation, branchEndPattern);
  // 
  //   int i = 0;
  //   if(local_debug > 0)
  //   {
  //     for(std::smatch::iterator it = matchLocation.begin(); it != matchLocation.end(); it++)
  //     {
  //       std::cout << i << ": " << *it << "\n";
  //       i++;
  //     }
  //   }
  // 
  //   std::string thisBranchString = matchLocation.str(0); //Need full string to get position
  //   std::size_t thisBranchEndBracket = queryStringCopy.find(thisBranchString);
  //   thisBranchEndBracket+=3;//Need location of the bracket
  //   char* endBracketPointer = &queryStringCopy[thisBranchEndBracket];
  //   std::size_t thisBranchStartBracket;
  //   unsigned int offset = queryStringCopy.size() - thisBranchEndBracket;//Offset for this position in reverse iterator
  //   int numOpenBrackets = 0;//Keeping track of [] in case of derivatives {IE [2P]}
  //   if(matchLocation.str(2) == "]")//found a branch
  //   {
  //     numOpenBrackets++;
  //     for(std::string::reverse_iterator rit = queryStringCopy.rbegin() + offset; rit != queryStringCopy.rend(); rit ++)
  //     {
  //       if(*rit == ']')
  //       {
  //         numOpenBrackets++;
  //       }
  //       else if(*rit == '[')
  //       {
  //         if(numOpenBrackets == 1)//This will close the bracket, and is the beginning of the branch
  //         {
  //           thisBranchStartBracket = std::distance(rit,queryStringCopy.rend()) - 1;//r(everse)end is the beginning
  //           numOpenBrackets--;
  //           break;
  //         }
  //         numOpenBrackets--;
  //       }
  //     }
  //     char* startBracketPointer = &queryStringCopy[thisBranchStartBracket];
  //     if(local_debug > 0)
  //     {
  //       std::cout << numOpenBrackets << " open brackets\n";
  //       std::cout << queryStringCopy << "\n";
  //       for(unsigned int i = 0; i < queryStringCopy.size(); i++)
  //       {
  //         if(i == thisBranchStartBracket)
  //         {
  //           std::cout << "[";
  //         }
  //         else if (i == thisBranchEndBracket)
  //         {
  //           std::cout << "]";
  //         }
  //         else
  //         {
  //           std::cout << " ";
  //         }
  //       }
  //       std::cout << '\n';
  //       for(unsigned int i = 0; i < queryStringCopy.size(); i++)
  //       {
  //         if(&queryStringCopy[i] == startBracketPointer)
  //         {
  //           std::cout << "[";
  //         }
  //         else if (&queryStringCopy[i] == endBracketPointer)
  //         {
  //           std::cout << "]";
  //         }
  //         else
  //         {
  //           std::cout << " ";
  //         }
  //       }
  //       std::cout << '\n';
  //       std::cout << &startBracketPointer << " " << &endBracketPointer << "\n";
  //     }
  //     std::string branchString = queryStringCopy.substr(thisBranchStartBracket + 1, thisBranchEndBracket - (thisBranchStartBracket + 1));
  //     queryStringCopy = queryStringCopy.substr(0,thisBranchStartBracket) + queryStringCopy.substr(thisBranchEndBracket+1, queryStringCopy.size());
  //     if(local_debug > 0)
  //     {
  //       std::cout << "Branch: " << branchString << " found\n";
  //       std::cout << queryStringCopy << "\n";
  //     }
  //   }
  //   else //No more branches
  //   {
  //     queryStringCopy = "";
  //   }
  // 
  // 
  // }
  // 
  
  return graph;
}

void MolecularModeling::Assembly::ConnectNodes(int start, int end, std::vector<parsedString>& parsedVector, GraphDS::Graph& graph)
{
  int local_debug = 1;
  if(local_debug > 0)
  {
    std::cout << "Connecting nodes\n";
    std::cout << start << ":" << end << "\n";
    std::cout << parsedVector.size()<< "\n";
  }
  
  if(end <= parsedVector.size())
  {
    for(unsigned int i=start; i<end; i++)
    {
      if(parsedVector[i].node != NULL)
      {
        graph.AddNewNode(parsedVector[i].node);
        if(i < parsedVector.size() - 2)
        {
          if(parsedVector[i+1].edge != NULL)
          {
            if(parsedVector[i+2].node != NULL)
            {
              parsedVector[i+1].edge->SetSourceNode(parsedVector[i].node);
              parsedVector[i+1].edge->SetDestinationNode(parsedVector[i+2].node);
              if(local_debug > 0)
              {
                std::cout << parsedVector[i].label << "{" << i <<"}" << " " << parsedVector[i+1].label << "{" << i + 1 <<"}" << " " <<  parsedVector[i + 2].label << "{" << i + 2 <<"}" << "\n";
              }
            }
            else if(parsedVector[i+2].label == "[")
            {//there is a branch.  Find the end of the branch and connect to the next node
              int numOpenBrackets = 1;
              bool onSameBranchPoint = true;
              int numBranchesAtThisNode = 1;
              int endBracketLocation;
              std::vector<int> startBranchLocations;
              startBranchLocations.push_back(i+2);
              std::vector<int> endBranchNodesLocations;
              while (onSameBranchPoint)
              {
                for(unsigned int j=i+3; j < end; j++)
                {//this will pass over both branched branches [[]] and multiple branches at the same node [][][]
                  //and hopefully point to the last bracket      ^                                             ^
                  if(parsedVector[j].label == "[")
                  {//branched branches
                    numOpenBrackets++;
                  }
                  if((parsedVector[j].label == "]") && (numOpenBrackets == 1))
                  {
                    if((j < parsedVector.size() - 1) && (parsedVector[j+1].label == "["))
                    {
                      numBranchesAtThisNode++;
                      startBranchLocations.push_back(j+2);
                      if(parsedVector[j-2].node != NULL)
                      {
                        endBranchNodesLocations.push_back(j-2);
                      }
                    }
                    else
                    {
                      if(parsedVector[j-2].node != NULL)
                      {
                        endBranchNodesLocations.push_back(j-2);
                      }
                      endBracketLocation = j;
                      onSameBranchPoint = false;
                      break;
                    }
                  }
                  else if((parsedVector[j].label == "]") && (numOpenBrackets != 1))
                  {
                    numOpenBrackets--;
                  }
                }
              }
              if((endBracketLocation < parsedVector.size() - 1) && (parsedVector[endBracketLocation + 1].node != NULL))
              {//Deal with the node before all this branching nonesense
                parsedVector[i+1].edge->SetSourceNode(parsedVector[i].node);
                parsedVector[i+1].edge->SetDestinationNode(parsedVector[endBracketLocation + 1].node);
                if(local_debug > 0)
                {
                  std::cout << parsedVector[i].label << "{" << i <<"}"<< " " <<  parsedVector[i+1].label << "{" << i +1 <<"}" << " " <<  parsedVector[endBracketLocation + 1].label << "{" << endBracketLocation + 1 <<"}" << "\n";
                }
              }
              for(int j = 0; j < endBranchNodesLocations.size(); j++)
              {//Attach the end of the branch to the node, recursively call this function to connect nodes in branch and deal with branched branches
                if((parsedVector[endBranchNodesLocations[j]+1].edge != NULL) && (parsedVector[endBracketLocation + 1].node != NULL))
                {
                  parsedVector[endBranchNodesLocations[j]+1].edge->SetSourceNode(parsedVector[j].node);
                  parsedVector[endBranchNodesLocations[j]+1].edge->SetDestinationNode(parsedVector[endBracketLocation + 1].node);
                  if(local_debug > 0)
                  {
                    std::cout << parsedVector[j].label << "{" << j <<"}" << " " <<  parsedVector[j+1].label<< "{" << j+1 <<"}" << " " <<  parsedVector[endBracketLocation + 1].label<< "{" << endBracketLocation + 1 <<"}" << "\n";
                  }
                }
                ConnectNodes(startBranchLocations[j], endBranchNodesLocations[j], parsedVector, graph);
              }
            }
          }
        }
      }
      else if(parsedVector[i].edge != NULL)
      {
        graph.AddEdge(parsedVector[i].edge);
      }
    }
  }
}

// To test new queries, you should go to your dev site's online virtuoso query.
// To find the IP address of your Virtuoso database, run docker inspect (YOUR_USERNAME)_gw_virt | grep IPAddress
// Then go to that IP address, followed by :8890/sparql.  For me right now, that's http://172.16.3.8:8890/sparql
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
//TODO make query work for any permutation of DGalpb1-4DGlcpNAcb1-2DManpa1-3[DGlcpNAcb1-2DManpa1-6]DManpb1-4DGlcpNAcb1-4[LFucpa1-6]DGlcpNAcb1-ASN

