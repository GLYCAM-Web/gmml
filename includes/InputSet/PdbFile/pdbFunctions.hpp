#ifndef INCLUDES_INPUTSET_PDBFILE_PDBFUNCTIONS_HPP
#define INCLUDES_INPUTSET_PDBFILE_PDBFUNCTIONS_HPP

#include <sstream> // stringstream

#include "includes/InputSet/PdbFile/pdbResidueId.hpp" // residueId

namespace pdb
{
//struct ResidueId
//{
//    std::string residueName_;
//    std::string sequenceNumber_;
//    std::string insertionCode_;
//    std::string chainId_;
//    bool operator==( ResidueId& rhs)
//                    { return (
//                            (residueName_ == rhs.residueName_) &&
//                            (sequenceNumber_ == rhs.sequenceNumber_) &&
//                            (insertionCode_ == rhs.insertionCode_) &&
//                            (chainId_ == rhs.chainId_));}
//    bool operator!=( ResidueId& rhs)
//                    { return (
//                            (residueName_ != rhs.residueName_) ||
//                            (sequenceNumber_ != rhs.sequenceNumber_) ||
//                            (insertionCode_ != rhs.insertionCode_) ||
//                            (chainId_ != rhs.chainId_));}
//};
//ResidueId extractResidueId(const std::string &line);
}
#endif /* INCLUDES_INPUTSET_PDBFILE_PDBFILE_HPP */
