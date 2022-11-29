#include "includes/InputSet/PdbFile/pdbFile.hpp"
#include "includes/InputSet/PdbFile/pdbModel.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp"
#include <fstream>      // std::ifstream
#include <algorithm>    // std::find

using pdb::PdbFile;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFile::PdbFile()
{
    //this->Initialize();
    inFilePath_ = "GMML-Generated";
}
PdbFile::PdbFile(const std::string &pdbFilePath) : inFilePath_(pdbFilePath)
{   //ToDo have a generic file opener in codeUtils
    std::ifstream pdbFileStream(pdbFilePath);
    if(pdbFileStream.fail())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Could not open this file: " + pdbFilePath);
        throw std::runtime_error("PdbFile constructor could not open this file: " + pdbFilePath);
    }
    gmml::log(__LINE__, __FILE__,  gmml::INF, "File opened: " + pdbFilePath + ". Ready to parse!");
    this->ParseInFileStream(pdbFileStream);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Finished parsing " + pdbFilePath);
}
void PdbFile::ParseInFileStream(std::ifstream& pdbFileStream)
{
    std::cout << "Parsing inputfile\n";
    for (std::string line; std::getline(pdbFileStream, line); )
    {
        //std::cout << "Parsing the line: " << line << "\n";
        codeUtils::ExpandLine(line, pdb::iPdbLineLength);
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        std::vector<std::string> coordSectionCards {"MODEL", "ATOM", "ANISOU", "TER", "HETATM", "CONECT"};
        std::vector<std::string> databaseCards {"DBREF", "DBREF1", "DBREF2"};
        if(std::find(coordSectionCards.begin(), coordSectionCards.end(), recordName) != coordSectionCards.end())
        {
            std::stringstream recordSection = this->ExtractHeterogenousRecordSection(pdbFileStream, line, coordSectionCards);
            this->addAssembly(std::make_unique<PdbModel>(recordSection)); // addAssembly is inherited
        }
        else if(recordName == "HEADER")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            headerRecord_ = HeaderRecord(recordSection);
        }
        else if(recordName == "TITLE")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            titleRecord_ = TitleRecord(recordSection);
        }
        else if(recordName == "AUTHOR")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            authorRecord_ = AuthorRecord(recordSection);
        }
        else if(recordName == "JRNL")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            journalRecord_ = JournalRecord(recordSection);
        }
        else if(recordName == "REMARK")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            remarkRecord_ = RemarkRecord(recordSection);
        }
        else if( (std::find(databaseCards.begin(), databaseCards.end(), recordName) != databaseCards.end()) )
        {
            std::stringstream databaseSection = this->ExtractHeterogenousRecordSection(pdbFileStream, line, databaseCards);
            while(getline(databaseSection, line))
            {
                databaseReferences_.emplace_back(line);
            }
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::INF, "PdbFile Constructor Complete Captain");
    return;
}
// Initializers used by constructors
// Should extract all lines that start with the strings in recordNames.
// Returns when it hits a line that does not start with one of those records.
std::stringstream PdbFile::ExtractHeterogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, const std::vector<std::string> recordNames)
{
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::stringstream recordSection;
    std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
    while(std::find(recordNames.begin(), recordNames.end(), recordName) != recordNames.end())
    {
        if(recordName != "ANISOU") // Do nothing for ANISOU
        {
            std::stringstream partialRecordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            recordSection << partialRecordSection.str();
            previousLinePosition = pdbFileStream.tellg(); // Save current line position.
        }
        if(!std::getline(pdbFileStream, line)) // If we hit the end
        {
            break; // Time to leave.
        }
        recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
    }
    pdbFileStream.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    //std::cout << "Returning this hetero section:\n" << heteroRecordSection.str() << "\nThe End of Hetero Section.\n";
    return recordSection;
}

// Goes through a section of the PDB file that contains the same header section. E.g. HEADER.
// If the header changes, it goes back to the previous line. I wanted the out while loop to trigger the new line. This means I don't have to check recordName between if statement and can have else if.
std::stringstream PdbFile::ExtractHomogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, std::string recordName)
{
    std::stringstream recordSection;
    codeUtils::ExpandLine(line, pdb::iPdbLineLength);
    recordSection << line << std::endl;
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::string previousName = recordName;
    while ( (std::getline(pdbFileStream, line)) )
    {
        codeUtils::ExpandLine(line, pdb::iPdbLineLength);
        recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        if(recordName != "ANISOU") // Do nothing for ANISOU
        {
            if (recordName == previousName)
            {
                recordSection << line << std::endl;
                previousName = recordName;
                previousLinePosition = pdbFileStream.tellg(); // Save current line position.
            }
            else
            {
                break;
            }
        } // Do nothing for ANISOU
    }
    pdbFileStream.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    //std::cout << "At end. Returning this record section:\n" << recordSection.str() << "\nEND RECORD SECTION\n";
    return recordSection;
}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
std::string PdbFile::GetUniprotIDs() const
{
    std::string UniprotIDs = "";
    for (auto &databaseReference : this->GetDatabaseReferences())
    {
        UniprotIDs += databaseReference.GetUniprotID();
    }
    return UniprotIDs;
}

const float& PdbFile::GetResolution() const
{
    return this->GetRemarkRecord().GetResolution();
}

const float& PdbFile::GetBFactor() const
{
    return this->GetRemarkRecord().GetBFactor();
}

pdb::PreprocessorInformation PdbFile::PreProcess(PreprocessorOptions inputOptions)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Preprocesssing has begun");
    pdb::PreprocessorInformation ppInfo;
    for(auto &cdsAssembly: this->getAssemblies()) // Now we do all, but maybe user can select at some point.
    {
        PdbModel* model = static_cast<PdbModel*>(cdsAssembly);
        model->preProcessCysResidues(ppInfo);
        model->preProcessHisResidues(ppInfo, inputOptions);
        model->preProcessChainTerminals(ppInfo, inputOptions);
        model->preProcessGapsUsingDistance(ppInfo, inputOptions);
        model->preProcessMissingUnrecognized(ppInfo);
        // ToDo What about setting charges doofus?
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Preprocessing completed");
    return ppInfo;
}

void PdbFile::Write(const std::string outName) const
{
    std::ofstream outFileStream;
    try
    {
        outFileStream.open(outName.c_str());
    }
    catch(...)
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Output file could not be created:\n" + outName);
        throw std::runtime_error("Output file could not be created:\n" + outName);
    }
    try
    {
        this->Write(outFileStream);
    }
    catch(...)
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Error when writing pdbFile class to file:\n" + outName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + outName);
    }
    return;
}

void PdbFile::Write(std::ostream& out) const
{
    this->GetHeaderRecord().Write(out);
    this->GetTitleRecord().Write(out);
    this->GetAuthorRecord().Write(out);
    this->GetJournalRecord().Write(out);
    this->GetRemarkRecord().Write(out);
    for (auto &dbref : this->GetDatabaseReferences())
    {
        dbref.Write(out);
    }
    for (auto &model : this->getAssemblies())
    {
        static_cast<PdbModel*>(model)->Write(out);
    }
    return;
}
