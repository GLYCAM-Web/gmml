// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfideresidue.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDisulfideResidue::PdbDisulfideResidue() : residue_name_(""), residue_chain_identifier_(' '), residue_insertion_code_(' '), symmetry_operator_(kNotSet) {}

PdbDisulfideResidue::PdbDisulfideResidue(const string &residue_name, char residue_chain_identifier, char residue_insertion_code, int symmetry_operator)
    : residue_name_(residue_name), residue_chain_identifier_(residue_chain_identifier), residue_insertion_code_(residue_insertion_code),
      symmetry_operator_(symmetry_operator) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbDisulfideResidue::GetResidueName()
{
    return residue_name_;
}

char PdbDisulfideResidue::GetResidueChainIdentifier()
{
    return residue_chain_identifier_;
}

char PdbDisulfideResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
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

void PdbDisulfideResidue::SetResidueChainIdentifier(char residue_chain_identifier)
{
    residue_chain_identifier_ = residue_chain_identifier;
}

void PdbDisulfideResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
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

