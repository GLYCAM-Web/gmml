#include "includes/InputSet/PdbFile/pdbResidueId.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/common.hpp"

using pdb::ResidueId;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
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

ResidueId::ResidueId(const std::string &name, const std::string &number, const std::string &insertionCode, const std::string &chainId) : residueName_(name), sequenceNumber_(number), insertionCode_(insertionCode), chainId_(chainId) {}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
const std::string ResidueId::getNumberAndInsertionCode() const
{
    return this->getNumber() + this->getInsertionCode();
}
//////////////////////////////////////////////////////////
//                       DISPLAY                        //
//////////////////////////////////////////////////////////
std::string ResidueId::print() const
{
    std::string formattedId;
    if (this->getName().empty())
    {
        formattedId += gmml::sNotSet;
    }
    else
    {
        formattedId += this->getName();
    }
    formattedId += "_";
    if (this->getNumberAndInsertionCode().empty())
    {
        formattedId += gmml::sNotSet;
    }
    else
    {
        formattedId += this->getNumberAndInsertionCode();
    }
    formattedId += "_";
    if (this->getChainId().empty())
    {
        formattedId += gmml::sNotSet;
    }
    else
    {
        formattedId += this->getChainId();
    }
    return formattedId;
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