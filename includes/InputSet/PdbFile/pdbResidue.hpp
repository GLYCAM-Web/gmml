#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP

#include <string>
#include <iostream>
#include <functional>
#include "includes/InputSet/PdbFile/atomRecord.hpp"

// Oliver Jan 2022
// This class just holds (non-owning) references to atomEntries so that I can more efficiently pass information
// out of the PdbFile class to other parts of GMML in a structure (residue) that is useful to the other parts.
// An equivalent Record to "residue" doesn't directly exist in the PDB format. Rather residues are made from the information in ATOM records. A group of ATOM records with the same "chainID" and "resSeq" (aka residue number) belong in the same "residue".
// This class is owned by CoordinateSection, which also owns the ATOM records. So the refs stored here should never be "dead".
// By ownership I mean responsible for creating and managing lifetime of.
namespace pdb
{
typedef std::vector<std::unique_ptr<AtomRecord>>::iterator AtomRecordIterator;
class PdbResidue : public abstrab::Labels
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
//    PdbResidue(AtomRecord* atomRecord);
//    PdbResidue(std::vector<AtomRecord*> atomRecords);
    PdbResidue(const std::string &line, const int& currentModelNumber);
    PdbResidue(const std::string residueName, const std::string atomName, GeometryTopology::Coordinate& atomCoord, const PdbResidue *referenceResidue);

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    std::string GetId() const;
    const std::string& GetChainId() const;
    const std::string& GetName() const;
    const std::string& GetRecordName() const;
    const std::string& GetInsertionCode() const;
    const int& GetSequenceNumber() const;
    std::vector<std::string> GetAtomNames() const;
    const std::string GetParmName() const;
    const int& GetModelNumber() const;
    AtomRecord* GetLastAtom() const;
    AtomRecord* GetFirstAtom() const;
    const std::string& GetLabel() const;
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //void AddAtom(AtomRecord* atomRecord);
    void CreateAtom(const std::string& line, const int& currentModelNumber);
    void CreateAtom(const std::string atomName, GeometryTopology::Coordinate& atomCoord);
    void SetName(const std::string name);
    inline void AddTerCard() {hasTerCard_ = true;}
    inline void RemoveTerCard() {hasTerCard_ = false;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    AtomRecord* FindAtom(const std::string& queryName) const;
    AtomRecord* FindAtom(const int& serialNumber) const;
    bool DeleteAtomRecord(AtomRecord* atom);
    AtomRecordIterator FindPositionOfAtom(AtomRecord* queryAtom);
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
    std::vector<std::unique_ptr<AtomRecord>> atomRecords_; // Residue does not own these. Owned by residues's owner.
    bool hasTerCard_;
    int modelNumber_;
};
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP
