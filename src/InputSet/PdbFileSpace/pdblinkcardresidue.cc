// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbLinkCardResidue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbLinkCardResidue::PdbLinkCardResidue() : atom_name_(""), alternate_location_indicator_(' '), residue_name_(""), residue_chain_id_(' '),
    residue_sequence_number_(gmml::iNotSet), residue_insertion_code_(' '), symmetry_operator_(gmml::iNotSet) {}

PdbLinkCardResidue::PdbLinkCardResidue(const std::string &atom_name, char alternate_location_indicator, const std::string &residue_name,
                               char residue_chain_id, int residue_sequence_number, char residue_insertion_code, int symmetry_operator)
    : atom_name_(atom_name), alternate_location_indicator_(alternate_location_indicator), residue_name_(residue_name),
      residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number), residue_insertion_code_(residue_insertion_code),
      symmetry_operator_(symmetry_operator) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbLinkCardResidue::GetAtomName()
{
    return atom_name_;
}

char PdbLinkCardResidue::GetAlternateLocationIndicator()
{
    return alternate_location_indicator_;
}

std::string PdbLinkCardResidue::GetResidueName()
{
    return residue_name_;
}

char PdbLinkCardResidue::GetResidueChainId()
{
    return residue_chain_id_;
}

int PdbLinkCardResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}

char PdbLinkCardResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}

int PdbLinkCardResidue::GetSymmetryOperator()
{
    return symmetry_operator_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbLinkCardResidue::SetAtomName(const std::string atom_name)
{
    atom_name_ = atom_name;
}

void PdbLinkCardResidue::SetAlternateLocationIndicator(char alternate_location_indicator)
{
    alternate_location_indicator_ = alternate_location_indicator;
}

void PdbLinkCardResidue::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}

void PdbLinkCardResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}

void PdbLinkCardResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}

void PdbLinkCardResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}

void PdbLinkCardResidue::SetSymmetryOperator(int symmetry_operator)
{
    symmetry_operator_ = symmetry_operator;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbLinkCardResidue::Print(std::ostream &out)
{
    out << "Atom Name: " << atom_name_
        << ", Alterante Location Identifier: " << alternate_location_indicator_
        << ", Residue Name: " << residue_name_
        << ", Residue Chain Identifier: " << residue_chain_id_
        << ", Residue Sequence Number: ";
    if(residue_sequence_number_ != gmml::iNotSet)
        out << residue_sequence_number_;
    else
        out << " ";
    out << ", Residue Insertion Code: " << residue_insertion_code_
        << ", Symmetry Operator: ";
    if(symmetry_operator_ != gmml::iNotSet)
        out << symmetry_operator_;
    else
        out << " ";
    out << std::endl << std::endl;
}
