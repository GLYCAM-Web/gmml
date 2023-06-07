#include "includes/Glycan/ontologyGmmlInterface.hpp"
#include "includes/Glycan/ontologyvocabulary.hpp"
#include "includes/utils.hpp"

void Ontology::PrintOntology(std::stringstream& ont_stream, pdb::PdbFile const &pdbFile)
{
    //Match formatting of Ontology
    std::stringstream uri;
    uri << Ontology::ONT_PREFIX;
    if(!pdbFile.GetHeaderRecord().GetIdentifierCode().empty())
    {
        uri << pdbFile.GetHeaderRecord().GetIdentifierCode();
    }
    else
    {
        std::string file = gmml::Split(pdbFile.GetInputFilePath().substr(pdbFile.GetInputFilePath().find_last_of('/') + 1), ".").at(0);
        // std::transform(file.begin(), file.end(),file.begin(), ::toupper);
        uri << file;
    }
    std::string uriStr = uri.str();
    std::transform(uriStr.begin(), uriStr.end(), uriStr.begin(), ::tolower);
    gmml::AddLiteral( uriStr, Ontology::TYPE, Ontology::PDB, ont_stream );
    //Return PDB_ID
    if(!pdbFile.GetHeaderRecord().GetIdentifierCode().empty())
        gmml::AddLiteral( uriStr, Ontology::id, pdbFile.GetHeaderRecord().GetIdentifierCode(), ont_stream );
    //Return Protein Acession Number
    //TODO add check to make sure that it is the Uniprot database reference
    if(!pdbFile.GetUniprotIDs().empty())
        gmml::AddLiteral( uriStr, Ontology::hasProteinID, pdbFile.GetUniprotIDs(), ont_stream );
    //Return Title
    if(!pdbFile.GetTitleRecord().GetTitle().empty())
        gmml::AddLiteral( uriStr, Ontology::hasTitle, pdbFile.GetTitleRecord().GetTitle(), ont_stream );
    //Return Authors
    if(!pdbFile.GetAuthorRecord().GetAuthor().empty())
        gmml::AddLiteral( uriStr, Ontology::hasAuthors, pdbFile.GetAuthorRecord().GetAuthor(), ont_stream );
    //Return Journal
    if(!pdbFile.GetJournalRecord().GetReference().empty())
    {
        //Return JOURNAL
        gmml::AddLiteral( uriStr, Ontology::hasJournal, pdbFile.GetJournalRecord().GetReference(), ont_stream );
        //Return DOI
        gmml::AddLiteral( uriStr, Ontology::hasDOI, pdbFile.GetJournalRecord().GetDOI(), ont_stream );
        //Return PMID
        gmml::AddLiteral( uriStr, Ontology::hasPMID, pdbFile.GetJournalRecord().GetPMID(), ont_stream );
    }
    //Return Remarks
    if(pdbFile.GetResolution() != gmml::dNotSet)
    {
        //Return Resolution
        gmml::AddDecimal( uriStr, Ontology::hasResolution, pdbFile.GetResolution(), ont_stream );
    }
    if(pdbFile.GetBFactor() != gmml::dNotSet)
    {
        //Return B Factor
        gmml::AddDecimal( uriStr, Ontology::hasBFactor, pdbFile.GetBFactor(), ont_stream );
    }
    return;
}
