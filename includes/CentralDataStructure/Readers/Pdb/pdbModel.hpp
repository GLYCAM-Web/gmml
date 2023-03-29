#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/conectRecord.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include <vector>
#include <iostream>
namespace pdb
{
class PdbAtom;
class PdbResidue;
class PdbChain;
class PdbModel : public cds::Assembly
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    PdbModel();
    PdbModel(std::stringstream& stream_block);
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void ChangeResidueName(const std::string& selector, const std::string& newName);
    std::string extractChainId(const std::string &line);
    std::stringstream extractSingleChainFromRecordSection(std::stringstream &stream_block, std::string line, const std::string& initialChainID);
    //Preprocessing functions
    void preProcessCysResidues(pdb::PreprocessorInformation &ppInfo);
    void preProcessHisResidues(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions);
    void preProcessChainTerminals(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions);
    void preProcessGapsUsingDistance(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions);
    void preProcessMissingUnrecognized(pdb::PreprocessorInformation &ppInfo);
    void bondAtomsByDistance();
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    void Print(std::ostream& out = std::cerr) const;
    void Write(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                       PRIVATE FUNCTIONS              //
    //////////////////////////////////////////////////////////
    void addConectRecord(const cds::Atom* atom1, const cds::Atom* atom2);
    inline const std::vector<ConectRecord>& GetConectRecords() const {return conectRecords_;}
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    std::vector<ConectRecord> conectRecords_;
};
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP
