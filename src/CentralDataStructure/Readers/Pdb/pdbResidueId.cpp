#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/constants.hpp"

using pdb::ResidueId;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
ResidueId::ResidueId(const std::string& line)
{ // atom number can overrun, shifting everything to the right.
    int shift            = pdb::checkShiftFromSerialNumberOverrun(line);
    alternativeLocation_ = codeUtils::RemoveWhiteSpace(line.substr(
        16 + shift, 1)); // I created this to use it only when reading, so I can discard the additional entries.
    residueName_         = codeUtils::RemoveWhiteSpace(line.substr(17 + shift, 3));
    chainId_             = codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
    int secondShift      = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
    sequenceNumber_      = codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift));
    // Insertion code gets shifted right by every overrun in residue number.
    insertionCode_       = codeUtils::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
}

ResidueId::ResidueId(const std::string name, const std::string number, const std::string insertionCode,
                     const std::string chainId)
    : residueName_(name), sequenceNumber_(number), insertionCode_(insertionCode), chainId_(chainId)
{}

ResidueId::ResidueId(std::vector<std::string> inputVector)
{
    //    std::cout << "Inputs:\n";
    //    for(auto & input : inputVector)
    //    {
    //        std::cout << input << std::endl;
    //    }
    if (inputVector.size() < 4)
    {
        throw std::runtime_error("ResidueId cannot be constructed from inputs");
    }
    if (inputVector.at(0) != "?")
    {
        residueName_ = inputVector.at(0);
    }
    if (inputVector.at(1) != "?")
    {
        sequenceNumber_ = inputVector.at(1);
    }
    if (inputVector.at(2) != "?")
    {
        insertionCode_ = inputVector.at(2);
    }
    if (inputVector.at(3) != "?")
    {
        chainId_ = inputVector.at(3);
    }
    // Just leave alternativeLocation as empty.
}

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
    std::string formattedId = "";
    if (this->getName().empty())
    {
        formattedId += constants::sNotSet;
    }
    else
    {
        formattedId += this->getName();
    }
    formattedId += "_";
    if (this->getNumber().empty())
    {
        formattedId += constants::sNotSet;
    }
    else
    {
        formattedId += this->getNumber();
    }
    formattedId += "_";
    if (this->getInsertionCode().empty())
    {
        formattedId += constants::sNotSet;
    }
    else
    {
        formattedId += this->getInsertionCode();
    }
    formattedId += "_";
    if (this->getChainId().empty())
    {
        formattedId += constants::sNotSet;
    }
    else
    {
        formattedId += this->getChainId();
    }
    return formattedId;
}

// pdb::ResidueId pdb::extractResidueId(const std::string &line)
//{
//     // Dealing with number overruns for serialNumber and residueNumber
//     ResidueId residueId;
//         int shift = pdb::checkShiftFromSerialNumberOverrun(line);
//     residueId.residueName_ = codeUtils::RemoveWhiteSpace(line.substr(17 + shift, 3));
//     residueId.chainId_ = codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
//     int secondShift = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
//     residueId.sequenceNumber_ = codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift));
//     // Insertion code gets shifted right by every overrun in residue number.
//     residueId.insertionCode_ = codeUtils::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
//     return residueId;
// }
