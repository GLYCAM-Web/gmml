#ifndef INCLUDES_INPUTSET_PDBFILE_PDBFILE_HPP
#define INCLUDES_INPUTSET_PDBFILE_PDBFILE_HPP

// ToDo split preprocessor into separate class that inherits from this one or is friends?.
// ToDo Get rid of coordinate section and just have this class hold the records
// ToDo Ownership Hierarchy of PdbFile->Models->Chains->Residues->AtomRecords? This would solve the TER problem when reading in a tleap generated file where there isn't a chain ID, but there is a TER card.
// ToDo Warning about gaps being 1 residue.
// ToDo make more direct queries here instead of giving out HeaderRecord etc.
// ToDo ACE/NME between residues with same number but an insertion code.
// ToDo expand tests to include amber input and multiple models.

#include <string>
#include <vector>
#include <fstream>      // std::ifstream

#include "includes/InputSet/PdbFile/headerRecord.hpp"
#include "includes/InputSet/PdbFile/databaseReferenceRecord.hpp"
#include "includes/InputSet/PdbFile/titleRecord.hpp"
#include "includes/InputSet/PdbFile/authorRecord.hpp"
#include "includes/InputSet/PdbFile/journalRecord.hpp"
#include "includes/InputSet/PdbFile/remarkRecord.hpp"
#include "includes/Resolver/NewPdbPreprocessor/pdbPreprocessorInputs.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/InputSet/PdbFile/pdbChain.hpp"
#include "includes/InputSet/PdbFile/pdbModel.hpp"
#include "includes/CentralDataStructure/cdsEnsemble.hpp"
//#include "includes/InputSet/PdbFile/conectRecord.hpp"

namespace pdb
{
const int iPdbLineLength = 80;
class PdbModel;
class PdbChain;
class PbdResidue;
class AtomRecord;
class PdbFile : public cds::cdsEnsemble<PdbModel, PdbChain, PdbResidue, AtomRecord>
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
    inline std::vector<const PdbModel*> getModels() const {return this->getAssemblies();} // renaming cds inherited getter for niceness.
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
    //////////////////////////////////////////////////////////
    //                        DISPLAY                       //
    //////////////////////////////////////////////////////////
    //void Print(std::ostream& out = std::cout) const;
    void Write(const std::string outName) const;
    void Write(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    //inline std::vector<PdbModel*>getModels() {return this->getAssemblies();}
    void ParseInFileStream(std::ifstream& pdbFileStream);
    std::stringstream ExtractHeterogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, const std::vector<std::string> recordNames);
    std::stringstream ExtractHomogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, std::string previousName);
    inline const std::vector<DatabaseReference>& GetDatabaseReferences() const {return databaseReferences_;}
//    inline const std::vector<ConectRecord>& GetConectRecords() const {return conectRecords_;}
//    inline PdbModel& GetCoordinateSection() {return coordinateSection_;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
//    void AddConnection(AtomRecord* atom1, AtomRecord* atom2);
    //void DeleteAtomRecord(AtomRecord* atom);
//    void ModifyNTerminal(const std::string& type, PdbResidue* nTerminalResidue);
//    void ModifyCTerminal(const std::string& type, PdbResidue* residue);
//    void InsertCap(const PdbResidue& residue, const std::string& type);
    //////////////////////////////////////////////////////////
    //                        ATTRIBUTES                    //
    //////////////////////////////////////////////////////////
    std::string inFilePath_;
    HeaderRecord headerRecord_;
    TitleRecord titleRecord_;
    AuthorRecord authorRecord_;
    JournalRecord journalRecord_;
    RemarkRecord remarkRecord_;
//    PdbModel coordinateSection_;
    std::vector<DatabaseReference> databaseReferences_;
};
}
#endif /* INCLUDES_INPUTSET_PDBFILE_PDBFILE_HPP */
