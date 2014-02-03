// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbLinkResidue::PdbLinkResidue() : atom_name_(""), alternate_location_indicator_(' '), residue_name_(""), residue_chain_identifier_(' '),
    residue_sequence_number_(kNotSet), residue_insertion_code_(' '), symmetry_operator_(kNotSet) {}

PdbLinkResidue::PdbLinkResidue(const string &atom_name, char alternate_location_indicator, const string &residue_name,
                               char residue_chain_identifier, int residue_sequence_number, char residue_insertion_code, int symmetry_operator)
    : atom_name_(atom_name), alternate_location_indicator_(alternate_location_indicator), residue_name_(residue_name),
      residue_chain_identifier_(residue_chain_identifier), residue_sequence_number_(residue_sequence_number), residue_insertion_code_(residue_insertion_code),
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

char PdbLinkResidue::GetResidueChainIdentifier()
{
    return residue_chain_identifier_;
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

void PdbLinkResidue::SetResidueChainIdentifier(char residue_chain_identifier)
{
    residue_chain_identifier_ = residue_chain_identifier;
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
    out << "Atom Name: " << atom_name_ << ", Alterante Location Identifier: " << alternate_location_indicator_ << ", Residue Name: " <<
           residue_name_ << ", Residue Chain Identifier: " << residue_chain_identifier_ << ", Residue Sequence Number: " <<
           residue_sequence_number_ << ", Residue Insertion Code: " << residue_insertion_code_ << ", Symmetry Operator: " <<
           symmetry_operator_ << endl << endl;
}
