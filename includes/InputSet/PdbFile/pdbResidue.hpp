#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP

#include <string>
#include <iostream>
#include <functional>

#include "includes/CentralDataStructure/cdsResidue.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"

namespace pdb
{
//typedef std::vector<std::unique_ptr<AtomRecord>>::iterator AtomRecordIterator;
class PdbResidue : public cds::cdsResidue<AtomRecord>
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
    void CreateAtomFromLine(const std::string& line, const int& currentModelNumber);
    //void CreateAtom(const std::string atomName, GeometryTopology::Coordinate& atomCoord);
    void SetName(const std::string name);
    inline void AddTerCard() {hasTerCard_ = true;}
    inline void RemoveTerCard() {hasTerCard_ = false;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
//    AtomRecord* FindAtom(const std::string& queryName) const;
//    AtomRecord* FindAtom(const int& serialNumber) const;
    bool DeleteAtomRecord(AtomRecord* atom);
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
    bool hasTerCard_;
    int modelNumber_;
};
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBRESIDUE_HPP
