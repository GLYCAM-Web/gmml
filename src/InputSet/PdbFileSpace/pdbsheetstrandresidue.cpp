// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbSheetStrandResidue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSheetStrandResidue::PdbSheetStrandResidue() : residue_name_(""), residue_chain_id_(' '), residue_sequence_number_(gmml::dNotSet), residue_insertion_code_(' ') {}

PdbSheetStrandResidue::PdbSheetStrandResidue(const std::string &residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code)
    : residue_name_(residue_name), residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number),
      residue_insertion_code_(residue_insertion_code) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbSheetStrandResidue::GetResidueName()
{
    return residue_name_;
}

char PdbSheetStrandResidue::GetResidueChainId()
{
    return residue_chain_id_;
}

int PdbSheetStrandResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}

char PdbSheetStrandResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbSheetStrandResidue::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}

void PdbSheetStrandResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}

void PdbSheetStrandResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}

void PdbSheetStrandResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSheetStrandResidue::Print(std::ostream &out)
{
    out << "Residue Name: " << residue_name_
        << ", Residue Cahin Identifier: " << residue_chain_id_
        << ", Residue Sequence Number: ";
    if( residue_sequence_number_ != gmml::iNotSet)
        out << residue_sequence_number_;
    else
        out << " ";
    out << ", Residue Insertion Code: " << residue_insertion_code_
        << std::endl;
}
