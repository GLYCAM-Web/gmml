#ifndef SRC_INPUTSET_PDBFILE_PDBFILE_HPP
#define SRC_INPUTSET_PDBFILE_PDBFILE_HPP

#include <string>
#include <vector>
#include <fstream>      // std::ifstream


#include "includes/InputSet/PdbFile/coordinateSection.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/InputSet/PdbFile/conectRecord.hpp"
#include "includes/InputSet/PdbFile/headerRecord.hpp"
#include "includes/InputSet/PdbFile/databaseReferenceRecord.hpp"
#include "includes/InputSet/PdbFile/titleRecord.hpp"
#include "includes/InputSet/PdbFile/authorRecord.hpp"
#include "includes/InputSet/PdbFile/journalRecord.hpp"
#include "includes/InputSet/PdbFile/remarkRecord.hpp"

namespace pdb
{
const int iPdbLineLength = 80;
class PdbFile
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    PdbFile();
    PdbFile(const std::string &pdbFilePath);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline const HeaderRecord& GetHeaderRecord() const {return headerRecord_;}
    inline const TitleRecord& GetTitleRecord() const {return titleRecord_;}
    inline const AuthorRecord& GetAuthorRecord() const {return authorRecord_;}
    inline const JournalRecord& GetJournalRecord() const {return journalRecord_;}
    inline const RemarkRecord& GetRemarkRecord() const {return remarkRecord_;}
    inline const std::string& GetPath() const {return inFilePath_;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    std::string GetUniprotIDs() const;
    const float& GetResolution() const;
    const float& GetBFactor() const;
private:
    void ParseInFileStream(std::ifstream& pdbFileStream);
    std::stringstream ExtractHeterogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, const std::vector<std::string> recordNames);
    std::stringstream ExtractHomogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, std::string previousName);
    std::string GetExpandedLine(std::ifstream& pdbFileStream);
    void ParseLine(std::string &line);
    inline const std::vector<DatabaseReference>& GetDatabaseReferences() const {return databaseReferences_;}
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
