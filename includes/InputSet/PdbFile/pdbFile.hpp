#ifndef SRC_INPUTSET_PDBFILE_PDBFILE_HPP
#define SRC_INPUTSET_PDBFILE_PDBFILE_HPP

#include <string>
#include <vector>
#include <fstream>      // std::ifstream

#include "includes/InputSet/PdbFile/coordinateSection.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/InputSet/PdbFile/headerRecord.hpp"
#include "includes/InputSet/PdbFile/databaseReferenceRecord.hpp"
#include "includes/InputSet/PdbFile/titleRecord.hpp"
#include "includes/InputSet/PdbFile/authorRecord.hpp"
#include "includes/InputSet/PdbFile/journalRecord.hpp"
#include "includes/InputSet/PdbFile/remarkRecord.hpp"
#include "includes/Resolver/NewPdbPreprocessor/pdbPreprocessorInputs.hpp"
#include "includes/InputSet/PdbFile/conectRecord.hpp"
#include "includes/InputSet/PdbFile/pdbResidue.hpp"

namespace pdb
{
const int iPdbLineLength = 80;
class PdbFile
{
   // friend class PdbPreprocessor;
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    PdbFile();
    PdbFile(const std::string &pdbFilePath);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    std::vector<PdbResidue> GetResiduesWithName(const std::string name) const;
    inline std::string GetInputFilePath() const {return inFilePath_;}
    // These should be private and whatever info they give out should be directly queryable here.
    inline const HeaderRecord& GetHeaderRecord() const {return headerRecord_;}
    inline const TitleRecord& GetTitleRecord() const {return titleRecord_;}
    inline const AuthorRecord& GetAuthorRecord() const {return authorRecord_;}
    inline const JournalRecord& GetJournalRecord() const {return journalRecord_;}
    inline const RemarkRecord& GetRemarkRecord() const {return remarkRecord_;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    std::string GetUniprotIDs() const;
    const float& GetResolution() const;
    const float& GetBFactor() const;
    pdb::PreprocessorInformation PreProcess(PreprocessorOptions options);
    void Print(std::ostream& out = std::cout) const;
private:
    void ParseInFileStream(std::ifstream& pdbFileStream);
    std::stringstream ExtractHeterogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, const std::vector<std::string> recordNames);
    std::stringstream ExtractHomogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, std::string previousName);
    inline const std::vector<DatabaseReference>& GetDatabaseReferences() const {return databaseReferences_;}
    inline CoordinateSection& GetCoordinateSection() {return coordinateSection_;}
    void AddConnection(AtomRecord* atom1, AtomRecord* atom2);
    void DeleteAtomRecord(AtomRecord* atom);
    void ModifyNTerminal(const std::string& type, PdbResidue* nTerminalResidue);
    void ModifyCTerminal(const std::string& type, PdbResidue* residue);
    void InsertCap(const PdbResidue& residue, const std::string& type);
    //////////////////////////////////////////////////////////
    //                        ATTRIBUTES                    //
    //////////////////////////////////////////////////////////
    std::string inFilePath_;
    HeaderRecord headerRecord_;
    TitleRecord titleRecord_;
    AuthorRecord authorRecord_;
    JournalRecord journalRecord_;
    RemarkRecord remarkRecord_;
    CoordinateSection coordinateSection_;
    std::vector<DatabaseReference> databaseReferences_;
    std::vector<ConectRecord> conectRecords_;
};
}
#endif /* SRC_INPUTSET_PDBFILE_PDBFILE_HPP */
