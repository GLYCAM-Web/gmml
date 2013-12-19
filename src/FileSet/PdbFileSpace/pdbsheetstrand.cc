// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbsheetstrand.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSheetStrand::PdbSheetStrand() : sense_(FIRST_STRAND), current_atom_(""), previous_atom_("") {}

PdbSheetStrand::PdbSheetStrand(const SheetStrandResidueVector strand_residues, PdbSheetStrandSense sense,
                               const string &current_atom, const string &previous_atom)
    : sense_(sense), current_atom_(current_atom), previous_atom_(previous_atom)
{
    strand_residues_.clear();
    for(PdbSheetStrand::SheetStrandResidueVector::const_iterator it = strand_residues.begin(); it != strand_residues.end(); it++)
    {
        strand_residues_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
PdbSheetStrand::SheetStrandResidueVector PdbSheetStrand::GetStrandResidues()
{
    return strand_residues_;
}

PdbSheetStrandSense PdbSheetStrand::GetSense()
{
    return sense_;
}

string PdbSheetStrand::GetCurrentAtom()
{
    return current_atom_;
}

string PdbSheetStrand::GetPreviousAtom()
{
    return previous_atom_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbSheetStrand::SetStrandResidues(const SheetStrandResidueVector strand_residues)
{
    strand_residues_.clear();
    for(PdbSheetStrand::SheetStrandResidueVector::const_iterator it = strand_residues.begin(); it != strand_residues.end(); it++)
    {
        strand_residues_.push_back(*it);
    }
}

void PdbSheetStrand::AddStrandResidue(PdbSheetStrandResidue *strand_residue)
{
    strand_residues_.push_back(strand_residue);
}

void PdbSheetStrand::SetSense(PdbSheetStrandSense sense)
{
    sense_ = sense;
}

void PdbSheetStrand::SetCurrentAtom(const string current_atom)
{
    current_atom_ = current_atom;
}

void PdbSheetStrand::SetPreviousAtom(const string previous_atom)
{
    previous_atom_ = previous_atom;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

