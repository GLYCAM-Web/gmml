#include <fstream>      // std::ifstream
#include <algorithm>    // std::find
#include "includes/InputSet/PdbFile/pdbFile.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/InputSet/PdbFile/databaseReferenceRecord.hpp"
#include "includes/common.hpp" // gmml::dSulfurCutoff

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
{
    std::ifstream pdbFileStream(pdbFilePath);
    if(pdbFileStream.fail())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Could not open this file: " + pdbFilePath);
        throw std::runtime_error("PdbFile constructor could not open this file: " + pdbFilePath);
    }
    gmml::log(__LINE__, __FILE__,  gmml::INF, "File opened: " + pdbFilePath + ". Ready to parse!");
    this->ParseInFileStream(pdbFileStream);
    coordinateSection_.Print();
    std::cout << "Finsihed printing Coordinate section" << std::endl;
}
void PdbFile::ParseInFileStream(std::ifstream& pdbFileStream)
{
    std::cout << "Parsing inputfile\n";
    for (std::string line; std::getline(pdbFileStream, line); )
    {
        //std::cout << "Parsing the line: " << line << "\n";
        codeUtils::ExpandLine(line, pdb::iPdbLineLength);
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        std::vector<std::string> cards {"MODEL", "ATOM", "ANISOU", "TER", "HETATM", "ENDMDL"};
        if(std::find(cards.begin(), cards.end(), recordName) != cards.end())
        {
            std::stringstream recordSection = this->ExtractHeterogenousRecordSection(pdbFileStream, line, cards);
            coordinateSection_ = CoordinateSection(recordSection);
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
    }
    return;
}
// Initializers used by constructors
std::stringstream PdbFile::ExtractHeterogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, const std::vector<std::string> recordNames)
{
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::stringstream recordSection;
    std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
    while(std::find(recordNames.begin(), recordNames.end(), recordName) != recordNames.end())
    {
        std::stringstream partialRecordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
        recordSection << partialRecordSection.str();
        previousLinePosition = pdbFileStream.tellg(); // Save current line position.
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

pdb::PreprocessorInformation PdbFile::PreProcess(PreprocessorOptions options)
{
    std::cout << "Preprocesssing has begun\n";
    PreprocessorInformation ppInfo;
    // CYS
    std::vector<pdb::PdbResidue> cysResidues = this->GetCoordinateSection().FindResidues("CYS");
    std::vector<pdb::PdbResidue> cyxResidues = this->GetCoordinateSection().FindResidues("CYX");
    cysResidues.insert( cysResidues.end(), cyxResidues.begin(), cyxResidues.end()); // concatenates the vectors.
    if (cysResidues.empty())
    {
        std::cout << "No CYS or CYX residue detected in this structure\n";
    }
    for (auto &cysRes1 : cysResidues)
    {
        AtomRecord* sg1 = cysRes1.FindAtom("SG");
        for (auto &cysRes2 : cysResidues)
        {
            AtomRecord* sg2 = cysRes2.FindAtom("SG");
            if ( (sg1 != nullptr) && (sg2 != nullptr) )
            {
                double distance = sg1->CalculateDistance(sg2);
                if (distance < gmml::dSulfurCutoff && distance > 0.001)
                {
                    cysRes1.SetName("CYX");
                    cysRes2.SetName("CYX");
                    //this->AddConnection(cysRes1.FindAtom("SG"), cysRes2.FindAtom("SG"));
                    std::cout << "Bonding " << cysRes1.GetId() << " and " << cysRes2.GetId() << " with distance " << distance << "\n";
                }
            }
        }
    }
    // HIS
    for (auto &hisRes : this->GetCoordinateSection().FindResidues("HIS"))
    {
        // HID residue
        if ( (hisRes.FindAtom("HE2") == nullptr) && (hisRes.FindAtom("HD1") != nullptr) )
        {
            hisRes.SetName("HID");
        }
        // HIP residue
        else if ( (hisRes.FindAtom("HE2") != nullptr) && (hisRes.FindAtom("HD1") != nullptr) )
        {
            hisRes.SetName("HIP");
        }
        // HIE is default
        else
        {
            hisRes.SetName("HIE");
        }
    }
    // Chain terminations
//    std::vector<pdb::PdbResidue> residues = this->GetCoordinateSection().GetChains();
    return ppInfo;
}
