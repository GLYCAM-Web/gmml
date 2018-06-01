// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbHelixResidue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHelixResidue::PdbHelixResidue() : residue_name_(""), residue_chain_id_(' '), residue_sequence_number_(gmml::dNotSet), residue_insertion_code_(' ' ) {}
PdbHelixResidue::PdbHelixResidue(const std::string &residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code)
    : residue_name_(residue_name), residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number),
      residue_insertion_code_(residue_insertion_code) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHelixResidue::GetResidueName()
{
    return residue_name_;
}

char PdbHelixResidue::GetResidueChainId()
{
    return residue_chain_id_;
}

int PdbHelixResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}

char PdbHelixResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHelixResidue::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}

void PdbHelixResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}

void PdbHelixResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}

void PdbHelixResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHelixResidue::Print(std::ostream &out)
{
    out << "Residue Name: " << residue_name_
        << ", Residue Chain ID: " << residue_chain_id_
        << ", Residue Sequence Number: ";
    if(residue_sequence_number_ != gmml::iNotSet)
        out << residue_sequence_number_;
    else
        out << " ";
    out << ", Residue Insertion Code: " << residue_insertion_code_
        << std::endl;
}
