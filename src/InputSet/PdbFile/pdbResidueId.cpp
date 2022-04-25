#include "includes/InputSet/PdbFile/pdbResidueId.hpp"
#include "includes/CodeUtils/strings.hpp"

using pdb::ResidueId;

ResidueId::ResidueId(const std::string &line)
{
    int shift = codeUtils::GetSizeOfIntInString(line.substr(12));
    residueName_ = codeUtils::RemoveWhiteSpace(line.substr(17 + shift, 3));
    chainId_ = codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
    int secondShift = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
    sequenceNumber_ = codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift));
    // Insertion code gets shifted right by every overrun in residue number.
    insertionCode_ = codeUtils::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
}

std::string ResidueId::getId() const
{
    return this->getName() + "_" + this->getNumberWithCode() + "_" + this->getChainId();
}

const std::string ResidueId::getNumberWithCode() const
{
    return this->getNumber() + this->getInsertionCode();
}

//pdb::ResidueId pdb::extractResidueId(const std::string &line)
//{
//    // Dealing with number overruns for serialNumber and residueNumber
//    ResidueId residueId;
//    int shift = codeUtils::GetSizeOfIntInString(line.substr(12));
//    residueId.residueName_ = codeUtils::RemoveWhiteSpace(line.substr(17 + shift, 3));
//    residueId.chainId_ = codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
//    int secondShift = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
//    residueId.sequenceNumber_ = codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift));
//    // Insertion code gets shifted right by every overrun in residue number.
//    residueId.insertionCode_ = codeUtils::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
//    return residueId;
//}
