// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbresiduemodification.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueModification::PdbResidueModification() : id_code_(""), residue_name_(""), chain_identifier_(' '), sequence_number_(kNotSet),
    insertion_code_(' '), standard_residue_name_(""), dscr_("") {}
PdbResidueModification::PdbResidueModification(const string &id_code, const string &residue_name, char chain_identifier, int sequence_number,
                                               char insertion_code, const string &standard_residue_name, const string &dscr) :
    id_code_(id_code), residue_name_(residue_name), chain_identifier_(chain_identifier), sequence_number_(sequence_number), insertion_code_(insertion_code),
    standard_residue_name_(standard_residue_name), dscr_(dscr) {}

PdbResidueModification::PdbResidueModification(stringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        id_code_ = line.substr(7,4);
        residue_name_ = line.substr(12,3);
        if(line.substr(16,1) == " ")
            chain_identifier_ = ' ';
        else
            chain_identifier_ = ConvertString<char>(line.substr(16,1));
        sequence_number_ = ConvertString<int>(line.substr(18,4));
        if(line.substr(22,1) == " ")
            insertion_code_ = ' ';
        else
            insertion_code_ = ConvertString<char>(line.substr(22,1));
        standard_residue_name_ = line.substr(24,3);
        dscr_ = line.substr(29,41);

        getline(stream_block, line);
        temp = line;
    }
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbResidueModification::GetIdCode()
{
    return id_code_;
}

std::string PdbResidueModification::GetResidueName()
{
    return residue_name_;
}

char PdbResidueModification::GetChainIdentifier()
{
    return chain_identifier_;
}

int PdbResidueModification::GetSequenceNumber()
{
    return sequence_number_;
}

char PdbResidueModification::GetInsertionCode()
{
    return insertion_code_;
}

std::string PdbResidueModification::GetStandardResidueName()
{
    return standard_residue_name_;
}

std::string PdbResidueModification::GetDscr()
{
    return dscr_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueModification::SetIdCode(const string id_code)
{
    id_code_ = id_code;
}

void PdbResidueModification::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}

void PdbResidueModification::SetChainIdentifier(char chain_identifier)
{
    chain_identifier_ = chain_identifier;
}

void PdbResidueModification::SetSequenceNumber(int sequence_number)
{
    sequence_number_ = sequence_number;
}

void PdbResidueModification::SetInsertionCode(char insertion_code)
{
    insertion_code_ = insertion_code;
}

void PdbResidueModification::SetStandardResidueName(const string standard_residue_name)
{
    standard_residue_name_ = standard_residue_name;
}

void PdbResidueModification::SetDscr(const string dscr)
{
    dscr_ = dscr;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueModification::Print(ostream &out)
{
    out << "ID Code: " << id_code_ << ", Residue Name: " << residue_name_ << ", Chain Identifier: " << chain_identifier_ <<
           ", Sequence Number: " << sequence_number_ << ", Insertion Code: " << insertion_code_ << ", Standard Residue Name: " << standard_residue_name_ <<
           ", Description: " << dscr_ << endl;
}
