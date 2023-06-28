#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS__PDB_PDBFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS__PDB_PDBFILE_HPP
// ToDo split preprocessor into separate function.
// ToDo Warning about gaps being 1 residue.
// ToDo make more direct queries here instead of giving out HeaderRecord etc.
// ToDo ACE/NME between residues with same number but an insertion code.
#include "includes/CentralDataStructure/ensemble.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/authorRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/databaseReferenceRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/headerRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/journalRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/remarkRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/titleRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include <string>
#include <vector>

namespace pdb
{
    const int iPdbLineLength = 80;
    class PdbModel;

    class PdbFile : public cds::Ensemble
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        PdbFile();
        PdbFile(const std::string& pdbFilePath);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline std::string GetInputFilePath() const
        {
            return inFilePath_;
        }

        // ToDo These should be private and whatever info they give out should be directly queryable here.
        inline const pdb::HeaderRecord& GetHeaderRecord() const
        {
            return headerRecord_;
        }

        inline const pdb::TitleRecord& GetTitleRecord() const
        {
            return titleRecord_;
        }

        inline const pdb::AuthorRecord& GetAuthorRecord() const
        {
            return authorRecord_;
        }

        inline const pdb::JournalRecord& GetJournalRecord() const
        {
            return journalRecord_;
        }

        inline const pdb::RemarkRecord& GetRemarkRecord() const
        {
            return remarkRecord_;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string GetUniprotIDs() const;
        const float& GetResolution() const;
        const float& GetBFactor() const;
        pdb::PreprocessorInformation PreProcess(PreprocessorOptions options);
        //////////////////////////////////////////////////////////
        //                        DISPLAY                       //
        //////////////////////////////////////////////////////////
        void Write(const std::string outName) const;
        void Write(std::ostream& stream) const;

      private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline const std::vector<pdb::DatabaseReference>& GetDatabaseReferences() const
        {
            return databaseReferences_;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void ParseInFileStream(std::ifstream& pdbFileStream);
        std::stringstream ExtractHeterogenousRecordSection(std::ifstream& pdbFileStream, std::string& line,
                                                           const std::vector<std::string> recordNames);
        std::stringstream ExtractHomogenousRecordSection(std::ifstream& pdbFileStream, std::string& line,
                                                         std::string previousName);
        //////////////////////////////////////////////////////////
        //                        ATTRIBUTES                    //
        //////////////////////////////////////////////////////////
        std::string inFilePath_;
        pdb::HeaderRecord headerRecord_; // SWIG wants the pdb::
        pdb::TitleRecord titleRecord_;
        pdb::AuthorRecord authorRecord_;
        pdb::JournalRecord journalRecord_;
        pdb::RemarkRecord remarkRecord_;
        std::vector<pdb::DatabaseReference> databaseReferences_;
    };
} // namespace pdb
#endif /* INCLUDES_CENTRALDATASTRUCTURE_READERS__PDB_PDBFILE_HPP */
