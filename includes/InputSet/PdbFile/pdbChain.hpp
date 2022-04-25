#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBCHAIN_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBCHAIN_HPP

#include <string>
#include <iostream>
#include <functional>

#include "includes/CentralDataStructure/cdsMolecule.hpp"

namespace pdb
{
class AtomRecord;
class PdbResidue;
class PdbChain : public cds::cdsMolecule<PdbResidue, AtomRecord>
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    PdbChain(std::stringstream &stream_block, const std::string& chainId);
//    PdbChain(PdbResidue pdbResidue);
//    PdbChain(std::vector<PdbResidue> pdbResidues);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    const std::string& GetChainId() const;
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void InsertCap(const PdbResidue& refResidue, const std::string& type);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
  //  void Print(std::ostream& out = std::cerr) const;
    void Write(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    //std::string extractResidueId(const std::string &line);
    std::stringstream extractSingleResidueFromRecordSection(std::stringstream &pdbFileStream, std::string line);

    //PdbResidue& GetFirstResidue() const;
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    //std::vector<PdbResidue> pdbResidues_;
    std::string chainId_;
};
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBCHAIN_HPP
