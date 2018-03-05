// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

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

PdbSheetStrand::PdbSheetStrand(string &line)
{
    int sense;
    // if(line.substr(38, 2).c_str() == "  ")
    //     sense = iNotSet;
    // else
        sense = ConvertString<int>(line.substr(38,2).c_str());
    switch(sense)
    {
        case -1:
            sense_ = ANTI_PARALLEL;
            break;
        case 0:
            sense_ = FIRST_STRAND;
            break;
        case 1:
            sense_ = PARALLEL;
            break;
        default:
            sense_ = UnknownStrand;
            break;
    }

    string temp0;
    char temp1, temp2;
    int temp3;
    temp0 = line.substr(17,3);
    Trim(temp0);
    if(line.substr(21,1).compare(" ") == 0)
        temp1 = ' ';
    else
        temp1 = ConvertString<char>(line.substr(21, 1));
    if(line.substr(26,1).compare(" ") == 0)
        temp2 = ' ';
    else
        temp2 = ConvertString<char>(line.substr(26, 1));
    if(line.substr(22, 4).compare("    ") == 0)
        temp3 = iNotSet;
    else
        temp3 = ConvertString<int>(line.substr(22, 4));
    PdbSheetStrandResidue* initial_residue = new PdbSheetStrandResidue(temp0, temp1, temp3, temp2);

    temp0 = line.substr(28, 3);
    Trim(temp0);
    if(line.substr(32,1).compare(" ") == 0)
        temp1 = ' ';
    else
        temp1 = ConvertString<char>(line.substr(32, 1));
    if(line.substr(37,1).compare(" ") == 0)
        temp2 = ' ';
    else
        temp2 = ConvertString<char>(line.substr(37, 1));
    if(line.substr(33, 4).compare("    ") == 0)
        temp3 = iNotSet;
    else
        temp3 = ConvertString<int>(line.substr(33, 4));
    PdbSheetStrandResidue* terminal_residue = new PdbSheetStrandResidue(temp0, temp1, temp3, temp2);
    strand_residues_.push_back(initial_residue);
    strand_residues_.push_back(terminal_residue);
    if(sense != 0)
    {
        temp0 = line.substr(45, 3);
        Trim(temp0);
        if(line.substr(49,1).compare(" ") == 0)
            temp1 = ' ';
        else
            temp1 = ConvertString<char>(line.substr(49, 1));
        if(line.substr(54,1).compare(" ") == 0)
            temp2 = ' ';
        else
            temp2 = ConvertString<char>(line.substr(54, 1));
        if(line.substr(50, 4).compare("    ") == 0)
            temp3 = iNotSet;
        else
            temp3 = ConvertString<int>(line.substr(50, 4));
        PdbSheetStrandResidue* current_residue = new PdbSheetStrandResidue(temp0, temp1, temp3, temp2);

        temp0 = line.substr(60, 3);
        Trim(temp0);
        if(line.substr(64,1).compare(" ") == 0)
            temp1 = ' ';
        else
            temp1 = ConvertString<char>(line.substr(64, 1));
        if(line.substr(69,1).compare(" ") == 0)
            temp2 = ' ';
        else
            temp2 = ConvertString<char>(line.substr(69, 1));
        if(line.substr(65, 4).compare("    ") == 0)
            temp3 = iNotSet;
        else
            temp3 = ConvertString<int>(line.substr(65, 4));
        PdbSheetStrandResidue* previous_residue = new PdbSheetStrandResidue(temp0, temp1, temp3, temp2);

        strand_residues_.push_back(current_residue);
        strand_residues_.push_back(previous_residue);
        current_atom_ = line.substr(41, 4);
        Trim(current_atom_);
        previous_atom_ = line.substr(56, 4);
        Trim(previous_atom_);
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
void PdbSheetStrand::Print(ostream &out)
{
    out << "Strand Residues: " << endl;
    for(PdbSheetStrand::SheetStrandResidueVector::iterator it = strand_residues_.begin(); it != strand_residues_.end(); it++)
    {
        (*it)->Print(out);
        out << endl;
    }
    out << "Sense: ";
    if(sense_ != UnknownStrand)
        out << sense_;
    else
        out << " ";
    out << ", Current Atom: " << current_atom_
        << ", Previous Atom: " << previous_atom_
        << endl << endl;
}
