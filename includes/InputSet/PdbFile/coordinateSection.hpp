#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP

#include <vector>
#include <iostream>
#include "includes/InputSet/PdbFile/pdbChain.hpp"
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"

namespace pdb
{
typedef std::vector<std::unique_ptr<PdbResidue>>::iterator PdbResidueIterator;
class CoordinateSection
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    CoordinateSection();
    CoordinateSection(std::stringstream& stream_block);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    // This next func is a stopgap until ParamManager comes into being.
    PdbResidue* CreateNewResidue(const std::string residueName, const std::string atomName, GeometryTopology::Coordinate& atomCoord, const PdbResidue& referenceResidue);
    std::vector<std::vector<pdb::PdbResidue*>> GetModels() const;
    std::vector<std::vector<pdb::PdbResidue*>> GetProteinChains(); // Exposed for PdbFile. Hmm maybe this shouldn't be a separate class.
    std::vector<pdb::PdbResidue*> GetResidues() const; // Exposed for PdbFile. Hmm maybe this shouldn't be a separate class.
    std::vector<pdb::PdbResidue*> FindResidues(const std::string selector); // Uh oh, private parts are exposed. Wee-ooo wee-ooo.
    void ChangeResidueName(const std::string& selector, const std::string& newName);
    AtomRecord* FindAtom(const int& serialNumber) const; // Conect records
  //  void DeleteAtomRecord(AtomRecord* atom);
//    AtomRecordIterator CreateNewAtomRecord(std::string name, GeometryTopology::Coordinate& coord, AtomRecord* sisterAtom);
//    AtomRecordIterator CreateNewAtomRecord(const std::string& atomName, const std::string& residueName, const int& residueSequenceNumber, const GeometryTopology::Coordinate& coord, const std::string& chainId, const int& modelNumber, AtomRecordIterator previousAtomPosition);
//    AtomRecordIterator FindPositionOfAtom(AtomRecord* queryAtom);
    PdbResidueIterator FindPositionOfResidue(const PdbResidue* queryResidue);
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
    //inline pdb::PdbResidue* GetCurrentResidue() {return &(residues_.front());}
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //inline pdb::Residue* CreateNewResidue(AtomRecord* atomRecord) {residues_.emplace_back(atomRecord);}
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    //std::vector<std::unique_ptr<AtomRecord>> atomRecords_;
    std::vector<std::unique_ptr<PdbResidue>> residues_;
    //std::vector<PdbResidue*> residues_;
};
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_COORDINATESECTION_HPP
