// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogen.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"
#include <iostream>

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogen::PdbHeterogen() : heterogen_id_(""), chain_identifier_(' '), sequence_number_(dNotSet), insertion_code_(' '),
    number_of_heterogen_atoms_(dNotSet), dscr_("") {}

PdbHeterogen::PdbHeterogen(const string &heterogen_id, char chain_identifier, int sequence_number,
                           char insertion_code, int number_of_heterogen_atoms, const string &dscr)
    : heterogen_id_(heterogen_id), chain_identifier_(chain_identifier), sequence_number_(sequence_number), insertion_code_(insertion_code),
      number_of_heterogen_atoms_(number_of_heterogen_atoms), dscr_(dscr) {}

PdbHeterogen::PdbHeterogen(stringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        heterogen_id_ = line.substr(7,3);
        Trim(heterogen_id_);
        if(line.substr(12,1) == " ")
        {
            chain_identifier_ = ' ';
        }
        else
        {
            chain_identifier_ = ConvertString<char>(line.substr(12,1));
        }
        if(line.substr(13, 4) == "    ")
            sequence_number_ = iNotSet;
        else
            sequence_number_ = ConvertString<int>(line.substr(13,4));
        if(line.substr(17,1) == " ")
        {
            insertion_code_ = ' ';
        }
        else
        {
            insertion_code_ = ConvertString<char>(line.substr(17,1));
        }
        if(line.substr(20, 5) == "     ")
            number_of_heterogen_atoms_ = iNotSet;
        else
            number_of_heterogen_atoms_ = ConvertString<int>(line.substr(20,5));
        dscr_ = line.substr(30,40);

        getline(stream_block, line);
        temp = line;
    }
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogen::GetHeterogenId()
{
    return heterogen_id_;
}

char PdbHeterogen::GetChainId()
{
    return chain_identifier_;
}

int PdbHeterogen::GetSequenceNumber()
{
    return sequence_number_;
}

char PdbHeterogen::GetInsertionCode()
{
    return insertion_code_;
}

int PdbHeterogen::GetNumberOfHeterogenAtoms()
{
    return number_of_heterogen_atoms_;
}

string PdbHeterogen::GetDscr()
{
    return dscr_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogen::SetHeterogenId(const string heterogen_id)
{
    heterogen_id_ = heterogen_id;
}

void PdbHeterogen::SetChainId(char chain_identifier)
{
    chain_identifier_ = chain_identifier;
}

void PdbHeterogen::SetSequenceNumber(int sequence_number)
{
    sequence_number_ = sequence_number;
}

void PdbHeterogen::SetInsertionCode(char insertion_code)
{
    insertion_code_ = insertion_code;
}

void PdbHeterogen::SetNumberOfHeterogenAtoms(int number_of_heterogen_atoms)
{
    number_of_heterogen_atoms_ = number_of_heterogen_atoms;
}

void PdbHeterogen::SetDscr(const string dscr)
{
    dscr_ = dscr;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogen::Print(ostream &out)
{
    out << "Heterogen ID: " << heterogen_id_
        << ", Chain Identifier: " << chain_identifier_
        << ", Sequence Number: ";
    if(sequence_number_ != iNotSet)
        out << sequence_number_;
    else
        out << " ";
    out << ", Insertion Code: " << insertion_code_
        << ", Number of Heterogen Atoms: ";
    if(number_of_heterogen_atoms_ != iNotSet)
        out << number_of_heterogen_atoms_;
    else
        out << " ";
    out << ", Description: " << dscr_ << endl;
}
