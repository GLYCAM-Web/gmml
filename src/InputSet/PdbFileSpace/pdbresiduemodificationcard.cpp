// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbResidueModificationCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueModificationCard::PdbResidueModificationCard() : id_code_(""), residue_name_(""), chain_id_(' '), sequence_number_(gmml::dNotSet),
    insertion_code_(' '), standard_residue_name_(""), dscr_("") {}
PdbResidueModificationCard::PdbResidueModificationCard(const std::string &id_code, const std::string &residue_name, char chain_id, int sequence_number,
                                               char insertion_code, const std::string &standard_residue_name, const std::string &dscr) :
    id_code_(id_code), residue_name_(residue_name), chain_id_(chain_id), sequence_number_(sequence_number), insertion_code_(insertion_code),
    standard_residue_name_(standard_residue_name), dscr_(dscr) {}

PdbResidueModificationCard::PdbResidueModificationCard(std::stringstream& stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        id_code_ = line.substr(7,4);
        gmml::Trim(id_code_);
        residue_name_ = line.substr(12,3);
        gmml::Trim(residue_name_);
        if(line.substr(16,1) == " ")
            chain_id_ = ' ';
        else
            chain_id_ = gmml::ConvertString<char>(line.substr(16,1));
        if(line.substr(18, 4) == "    ")
            sequence_number_ = gmml::iNotSet;
        else
            sequence_number_ = gmml::ConvertString<int>(line.substr(18,4));
        if(line.substr(22,1) == " ")
            insertion_code_ = ' ';
        else
            insertion_code_ = gmml::ConvertString<char>(line.substr(22,1));
        standard_residue_name_ = line.substr(24,3);
        gmml::Trim(standard_residue_name_);
        dscr_ = line.substr(29,41);

        getline(stream_block, line);
        temp = line;
    }
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbResidueModificationCard::GetIdCode()
{
    return id_code_;
}

std::string PdbResidueModificationCard::GetResidueName()
{
    return residue_name_;
}

char PdbResidueModificationCard::GetChainId()
{
    return chain_id_;
}

int PdbResidueModificationCard::GetSequenceNumber()
{
    return sequence_number_;
}

char PdbResidueModificationCard::GetInsertionCode()
{
    return insertion_code_;
}

std::string PdbResidueModificationCard::GetStandardResidueName()
{
    return standard_residue_name_;
}

std::string PdbResidueModificationCard::GetDscr()
{
    return dscr_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueModificationCard::SetIdCode(const std::string id_code)
{
    id_code_ = id_code;
}

void PdbResidueModificationCard::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}

void PdbResidueModificationCard::SetChainId(char chain_id)
{
    chain_id_ = chain_id;
}

void PdbResidueModificationCard::SetSequenceNumber(int sequence_number)
{
    sequence_number_ = sequence_number;
}

void PdbResidueModificationCard::SetInsertionCode(char insertion_code)
{
    insertion_code_ = insertion_code;
}

void PdbResidueModificationCard::SetStandardResidueName(const std::string standard_residue_name)
{
    standard_residue_name_ = standard_residue_name;
}

void PdbResidueModificationCard::SetDscr(const std::string dscr)
{
    dscr_ = dscr;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueModificationCard::Print(std::ostream &out)
{
    out << "ID Code: " << id_code_
        << ", Residue Name: " << residue_name_
        << ", Chain Identifier: " << chain_id_
        << ", Sequence Number: ";
    if(sequence_number_ != gmml::iNotSet)
        out << sequence_number_;
    else
        out << " ";
    out << ", Insertion Code: " << insertion_code_
        << ", Standard Residue Name: " << standard_residue_name_
        << ", Description: " << dscr_ << std::endl;
}
