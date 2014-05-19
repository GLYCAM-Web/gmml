// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfideresidue.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDisulfideResidue::PdbDisulfideResidue() : residue_name_(""), residue_chain_id_(' '), residue_sequence_number_(dNotSet), residue_insertion_code_(' '), symmetry_operator_(dNotSet) {}

PdbDisulfideResidue::PdbDisulfideResidue(const string &residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code, int symmetry_operator)
    : residue_name_(residue_name), residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number), residue_insertion_code_(residue_insertion_code),
      symmetry_operator_(symmetry_operator) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbDisulfideResidue::GetResidueName()
{
    return residue_name_;
}

char PdbDisulfideResidue::GetResidueChainId()
{
    return residue_chain_id_;
}

char PdbDisulfideResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}

int PdbDisulfideResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}

int PdbDisulfideResidue::GetSymmetryOperator()
{
    return symmetry_operator_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbDisulfideResidue::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}

void PdbDisulfideResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}

void PdbDisulfideResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}

void PdbDisulfideResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}

void PdbDisulfideResidue::SetSymmetryOperator(int symmetry_operator)
{
    symmetry_operator_ = symmetry_operator;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbDisulfideResidue::Print(ostream &out)
{
    out << "Residue Name: " << residue_name_
        << ", Residue Chain Identifier: " << residue_chain_id_
        << ", Residue Insertion Code: " << residue_insertion_code_
        << ", Residue Sequence Number: ";
    if(residue_sequence_number_ != iNotSet)
        out << residue_sequence_number_;
    else
        out << " ";
    out << ", Symmetry Operator: ";
    if(symmetry_operator_ != iNotSet)
        out << symmetry_operator_;
    else
        out << " ";
    out << endl;
}
