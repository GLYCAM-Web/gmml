// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueModificationCard::PdbResidueModificationCard() : id_code_(""), residue_name_(""), chain_id_(' '), sequence_number_(dNotSet),
    insertion_code_(' '), standard_residue_name_(""), dscr_("") {}
PdbResidueModificationCard::PdbResidueModificationCard(const string &id_code, const string &residue_name, char chain_id, int sequence_number,
                                               char insertion_code, const string &standard_residue_name, const string &dscr) :
    id_code_(id_code), residue_name_(residue_name), chain_id_(chain_id), sequence_number_(sequence_number), insertion_code_(insertion_code),
    standard_residue_name_(standard_residue_name), dscr_(dscr) {}

PdbResidueModificationCard::PdbResidueModificationCard(stringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        id_code_ = line.substr(7,4);
        Trim(id_code_);
        residue_name_ = line.substr(12,3);
        Trim(residue_name_);
        if(line.substr(16,1) == " ")
            chain_id_ = ' ';
        else
            chain_id_ = ConvertString<char>(line.substr(16,1));
        if(line.substr(18, 4) == "    ")
            sequence_number_ = iNotSet;
        else
            sequence_number_ = ConvertString<int>(line.substr(18,4));
        if(line.substr(22,1) == " ")
            insertion_code_ = ' ';
        else
            insertion_code_ = ConvertString<char>(line.substr(22,1));
        standard_residue_name_ = line.substr(24,3);
        Trim(standard_residue_name_);
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
void PdbResidueModificationCard::SetIdCode(const string id_code)
{
    id_code_ = id_code;
}

void PdbResidueModificationCard::SetResidueName(const string residue_name)
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

void PdbResidueModificationCard::SetStandardResidueName(const string standard_residue_name)
{
    standard_residue_name_ = standard_residue_name;
}

void PdbResidueModificationCard::SetDscr(const string dscr)
{
    dscr_ = dscr;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueModificationCard::Print(ostream &out)
{
    out << "ID Code: " << id_code_
        << ", Residue Name: " << residue_name_
        << ", Chain Identifier: " << chain_id_
        << ", Sequence Number: ";
    if(sequence_number_ != iNotSet)
        out << sequence_number_;
    else
        out << " ";
    out << ", Insertion Code: " << insertion_code_
        << ", Standard Residue Name: " << standard_residue_name_
        << ", Description: " << dscr_ << endl;
}
