// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbLinkResidue::PdbLinkResidue() : atom_name_(""), alternate_location_indicator_(' '), residue_name_(""), residue_chain_id_(' '),
    residue_sequence_number_(iNotSet), residue_insertion_code_(' '), symmetry_operator_(iNotSet) {}

PdbLinkResidue::PdbLinkResidue(const string &atom_name, char alternate_location_indicator, const string &residue_name,
                               char residue_chain_id, int residue_sequence_number, char residue_insertion_code, int symmetry_operator)
    : atom_name_(atom_name), alternate_location_indicator_(alternate_location_indicator), residue_name_(residue_name),
      residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number), residue_insertion_code_(residue_insertion_code),
      symmetry_operator_(symmetry_operator) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbLinkResidue::GetAtomName()
{
    return atom_name_;
}

char PdbLinkResidue::GetAlternateLocationIndicator()
{
    return alternate_location_indicator_;
}

string PdbLinkResidue::GetResidueName()
{
    return residue_name_;
}

char PdbLinkResidue::GetResidueChainId()
{
    return residue_chain_id_;
}

int PdbLinkResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}

char PdbLinkResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}

int PdbLinkResidue::GetSymmetryOperator()
{
    return symmetry_operator_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbLinkResidue::SetAtomName(const string atom_name)
{
    atom_name_ = atom_name;
}

void PdbLinkResidue::SetAlternateLocationIndicator(char alternate_location_indicator)
{
    alternate_location_indicator_ = alternate_location_indicator;
}

void PdbLinkResidue::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}

void PdbLinkResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}

void PdbLinkResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}

void PdbLinkResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}

void PdbLinkResidue::SetSymmetryOperator(int symmetry_operator)
{
    symmetry_operator_ = symmetry_operator;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbLinkResidue::Print(ostream &out)
{
    out << "Atom Name: " << atom_name_
        << ", Alterante Location Identifier: " << alternate_location_indicator_
        << ", Residue Name: " << residue_name_
        << ", Residue Chain Identifier: " << residue_chain_id_
        << ", Residue Sequence Number: ";
    if(residue_sequence_number_ != iNotSet)
        out << residue_sequence_number_;
    else
        out << " ";
    out << ", Residue Insertion Code: " << residue_insertion_code_
        << ", Symmetry Operator: ";
    if(symmetry_operator_ != iNotSet)
        out << symmetry_operator_;
    else
        out << " ";
    out << endl << endl;
}
