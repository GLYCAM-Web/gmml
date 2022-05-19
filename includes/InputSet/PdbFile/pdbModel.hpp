#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP

#include <vector>
#include <iostream>

#include "includes/CentralDataStructure/cdsAssembly.hpp"
#include "includes/InputSet/PdbFile/conectRecord.hpp"
#include "includes/Resolver/NewPdbPreprocessor/pdbPreprocessorInputs.hpp"
#include "includes/InputSet/PdbFile/pdbChain.hpp"

namespace pdb
{
class pdbAtom;
class PdbResidue;
class PdbChain;
class PdbModel : public cds::cdsAssembly<PdbChain, PdbResidue, pdbAtom>
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    PdbModel();
    PdbModel(std::stringstream& stream_block);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void ChangeResidueName(const std::string& selector, const std::string& newName);
//    const pdbAtom* FindAtom(const int& serialNumber) const; // Conect records
    std::string extractChainId(const std::string &line);
    std::stringstream extractSingleChainFromRecordSection(std::stringstream &stream_block, std::string line, const std::string& initialChainID);
    //Preprocessing functions
    void preProcessCysResidues(pdb::PreprocessorInformation &ppInfo);
    void preProcessHisResidues(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions);
    void preProcessChainTerminals(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions);
    void preProcessGaps(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions);
    void preProcessMissingUnrecognized(pdb::PreprocessorInformation &ppInfo);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    void Print(std::ostream& out = std::cerr) const;
    void Write(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void addConectRecord(const pdbAtom* atom1, const pdbAtom* atom2);
    inline const std::vector<ConectRecord>& GetConectRecords() const {return conectRecords_;}
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    std::vector<ConectRecord> conectRecords_;
};
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP
