#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP

#include <string>
#include <iostream>
#include <functional>

#include "includes/InputSet/PdbFile/pdbResidueId.hpp"
#include "includes/CentralDataStructure/cdsResidue.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"

namespace pdb
{
//typedef std::vector<std::unique_ptr<AtomRecord>>::iterator AtomRecordIterator;
class PdbResidue : public cds::cdsResidue<AtomRecord> , public pdb::ResidueId
{
public:
    using ResidueId::getName; // I want ResidueId's implementation of these methods
    using ResidueId::getNumber;
    using ResidueId::setName;
    using ResidueId::setNumber;
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
//    PdbResidue(AtomRecord* atomRecord);
//    PdbResidue(std::vector<AtomRecord*> atomRecords);
    PdbResidue(std::stringstream &singleResidueSecion, std::string firstLine);
    PdbResidue(const std::string residueName, const PdbResidue *referenceResidue);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    //std::string GetId() const;
    //inline const std::string& getChainId() const {return chainId_;}
    //inline const std::string& getInsertionCode() const {return insertionCode_;}
    //const std::string& GetName() const;
    const std::string& GetRecordName() const;
    std::vector<std::string> GetAtomNames() const;
    const std::string GetParmName() const;
//    AtomRecord* GetLastAtom() const;
//    AtomRecord* GetFirstAtom() const;
    const std::string& GetLabel() const;
    //const std::string& GetId() const;
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //inline void setChainId(const std::string& s) {chainId_ = s;}
    //inline void setInsertionCode(const std::string& s) {insertionCode_ = s;}
    inline void AddTerCard() {hasTerCard_ = true;}
    inline void RemoveTerCard() {hasTerCard_ = false;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    //void CreateAtomFromLine(const std::string& line, const int& currentModelNumber);
    void modifyNTerminal(const std::string& type);
    void modifyCTerminal(const std::string& type);
//    AtomRecord* FindAtom(const std::string& queryName) const;
//    AtomRecord* FindAtom(const int& serialNumber) const;
   // bool DeleteAtomRecord(AtomRecord* atom);
    //AtomRecordIterator FindPositionOfAtom(AtomRecord* queryAtom);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    void Print(std::ostream& out = std::cerr) const;
    void Write(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    //std::vector<AtomRecord*> atomRecords_; // Residue does not own these. Owned by residues's owner.
    //std::vector<std::unique_ptr<AtomRecord>> atomRecords_; // Residue does not own these. Owned by residues's owner.
//    std::string chainId_ = "";
//    std::string insertionCode_ = "";
    bool hasTerCard_ = false;
};
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP
