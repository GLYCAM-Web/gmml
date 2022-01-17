#ifndef SRC_INPUTSET_PDBFILE_PDBFILE_HPP
#define SRC_INPUTSET_PDBFILE_PDBFILE_HPP

#include <string>
#include <vector>

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
class PdbFile
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    PdbFile(const std::string &pdbFile);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline const HeaderRecord& GetHeaderRecord() const {return &headerRecord_;}
    inline const TitleRecord& GetTitleRecord() const {return &titleRecord_;}
    inline const AuthorRecord& GetAuthorRecord() const {return &authorRecord_;}
    inline const JournalRecord& GetJournalRecord() const {return &journalRecord_;}
    inline const RemarkRecord& GetRemarkRecord() const {return &remarkRecord_;}
    inline const std::string& GetPath() const {return &inFilePath_;}
    std::string GetUniprotIDs() const;
private:
    inline const std::vector<DatabaseReference> GetDatabaseReferences() const {return &databaseReferences_;}
    //////////////////////////////////////////////////////////
    //                        ATTRIBUTES                    //
    //////////////////////////////////////////////////////////
    std::string inFilePath_;
    std::ifstream inFileStream_;
    HeaderRecord headerRecord_;
    TitleRecord titleRecord_;
    AuthorRecord authorRecord_;
    JournalRecord journalRecord_;
    RemarkRecord remarkRecord_;
    std::vector<DatabaseReference> databaseReferences_;
    std::vector<AtomRecord> atomEntries_;
    std::vector<AtomRecord> hetatmEntries_;
    std::vector<ConectRecord> conectRecords_;
};
}
#endif /* SRC_INPUTSET_PDBFILE_PDBFILE_HPP */
