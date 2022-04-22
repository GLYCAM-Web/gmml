#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP

#include <vector>
#include <iostream>
//#include "includes/InputSet/PdbFile/pdbChain.hpp"
//#include "includes/InputSet/PdbFile/pdbResidue.hpp"
//#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/CentralDataStructure/cdsAssembly.hpp"
#include "includes/InputSet/PdbFile/conectRecord.hpp"
#include "includes/Resolver/NewPdbPreprocessor/pdbPreprocessorInputs.hpp"

namespace pdb
{
class AtomRecord;
class PdbResidue;
class PdbChain;
//typedef std::vector<std::unique_ptr<PdbResidue>>::iterator PdbResidueIterator;
class PdbModel : public cds::cdsAssembly<PdbChain, PdbResidue, AtomRecord>
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
    inline int getModelNumber() const {return modelNumber_;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    // This next func is a stopgap until ParamManager comes into being.
//    PdbResidue* CreateNewResidue(const std::string residueName, const std::string atomName, GeometryTopology::Coordinate& atomCoord, const PdbResidue& referenceResidue);
    //std::vector<std::vector<pdb::PdbResidue*>> GetModels() const;
    std::vector<pdb::PdbChain*> GetProteinChains(); // Exposed for PdbFile. Hmm maybe this shouldn't be a separate class.
    //std::vector<pdb::PdbResidue*> GetResidues() const; // Exposed for PdbFile. Hmm maybe this shouldn't be a separate class.
    std::vector<pdb::PdbResidue*> FindResidues(const std::string selector); // Uh oh, private parts are exposed. Wee-ooo wee-ooo.
    void ChangeResidueName(const std::string& selector, const std::string& newName);
    const AtomRecord* FindAtom(const int& serialNumber) const; // Conect records
  //  void DeleteAtomRecord(AtomRecord* atom);
//    AtomRecordIterator CreateNewAtomRecord(std::string name, GeometryTopology::Coordinate& coord, AtomRecord* sisterAtom);
//    AtomRecordIterator CreateNewAtomRecord(const std::string& atomName, const std::string& residueName, const int& residueSequenceNumber, const GeometryTopology::Coordinate& coord, const std::string& chainId, const int& modelNumber, AtomRecordIterator previousAtomPosition);
//    AtomRecordIterator FindPositionOfAtom(AtomRecord* queryAtom);
    //PdbResidueIterator FindPositionOfResidue(const PdbResidue* queryResidue);
    std::string extractChainId(const std::string &line);
    std::stringstream extractSingleChainFromRecordSection(std::stringstream &pdbFileStream, std::string line, const std::string& initialChainID);
    void preProcessCysResidues(pdb::PreprocessorInformation &ppInfo, const PreprocessorOptions &inputOptions);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    void Print(std::ostream& out = std::cerr) const;
    void Write(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    std::string PeekAtResidueId(const std::string &line);
    void addConectRecord(const AtomRecord* atom1, const AtomRecord* atom2);
    inline const std::vector<ConectRecord>& GetConectRecords() const {return conectRecords_;}
    //inline pdb::PdbResidue* GetCurrentResidue() {return &(residues_.front());}
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //inline pdb::Residue* CreateNewResidue(AtomRecord* atomRecord) {residues_.emplace_back(atomRecord);}
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    //std::vector<std::unique_ptr<AtomRecord>> atomRecords_;
   // std::vector<std::unique_ptr<PdbResidue>> residues_;
    //std::vector<PdbResidue*> residues_;
    int modelNumber_ = 1;
    std::vector<ConectRecord> conectRecords_;
};
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBMODEL_HPP
