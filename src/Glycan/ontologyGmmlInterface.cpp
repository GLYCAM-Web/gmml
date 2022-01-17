#include "includes/Glycan/ontologyGmmlInterface.hpp"
#include "includes/Glycan/ontologyvocabulary.hpp"
#include "includes/utils.hpp"

void Ontology::PrintOntology(std::stringstream& ont_stream, pdb::PdbFile const &pdbFile)
{
    //Match formatting of Ontology
    std::stringstream uri;
    uri << Ontology::ONT_PREFIX;
    if(pdbFile.GetHeaderRecord() != NULL)
    {
        uri << pdbFile.GetHeaderRecord().GetIdentifierCode();
    }
    else
    {
        std::string file = gmml::Split(pdbFile.GetPath().substr(pdbFile.GetPath().find_last_of('/') + 1), ".").at(0);
        // std::transform(file.begin(), file.end(),file.begin(), ::toupper);
        uri << file;
    }
    std::string uriStr = uri.str();
    std::transform(uriStr.begin(), uriStr.end(), uriStr.begin(), ::tolower);

    gmml::AddLiteral( uriStr, Ontology::TYPE, Ontology::PDB, ont_stream );

    //Return PDB_ID
    if(pdbFile.GetHeaderRecord() != NULL)
        gmml::AddLiteral( uriStr, Ontology::id, pdbFile.GetHeaderRecord().GetIdentifierCode(), ont_stream );

    //Return Protein Acession Number
    //TODO add check to make sure that it is the Uniprot database reference
    if(pdbFile.GetUniprotIDs() != "")
        gmml::AddLiteral( uriStr, Ontology::hasProteinID, pdbFile.GetUniprotIDs(), ont_stream );

    //Return Title
    if(pdbFile.GetTitleRecord() != NULL)
        gmml::AddLiteral( uriStr, Ontology::hasTitle, pdbFile.GetTitleRecord().GetTitle(), ont_stream );

    //Return Authors
    if(pdbFile.GetAuthorRecord() != NULL)
        gmml::AddLiteral( uriStr, Ontology::hasAuthors, pdbFile.GetAuthorRecord().GetAuthor(), ont_stream );


    if(pdbFile.GetJournalRecord() != NULL)
    {
        //Return JOURNAL
        gmml::AddLiteral( uriStr, Ontology::hasJournal, pdbFile.GetJournalRecord().GetReference(), ont_stream );

        //Return DOI
        gmml::AddLiteral( uriStr, Ontology::hasDOI, pdbFile.GetJournalRecord().GetDOI(), ont_stream );

        //Return PMID
        gmml::AddLiteral( uriStr, Ontology::hasPMID, pdbFile.GetJournalRecord().GetPMID(), ont_stream );
    }
    if(pdbFile.GetRemarkRecord() != NULL)
    {
        //Return Resolution
        gmml::AddDecimal( uriStr, Ontology::hasResolution, pdbFile.GetRemarkRecord().GetResolution(), ont_stream );

        //Return B Factor
        gmml::AddDecimal( uriStr, Ontology::hasBFactor, pdbFile.GetRemarkRecord().GetBFactor(), ont_stream );
    }
    return;
}
