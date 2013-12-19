// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbhelix.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHelix::PdbHelix() : helix_id_(""), helix_serial_number_(kNotSet), helix_class_(POLYPROLINE), comment_(""), helix_length_(kNotSet) {}

PdbHelix::PdbHelix(const string &helix_id, int helix_serial_number, HelixResidueVector helix_residues,
                   PdbHelixClass helix_class, const string &comment, double helix_length)
    : helix_id_(helix_id), helix_serial_number_(helix_serial_number), helix_class_(helix_class), comment_(comment), helix_length_(helix_length)
{
    helix_residues_.clear();
    for(PdbHelix::HelixResidueVector::const_iterator it = helix_residues.begin(); it != helix_residues.end(); it++)
    {
        helix_residues_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHelix::GetHelixId()
{
    return helix_id_;
}

int PdbHelix::GetHelixSerialNumber()
{
    return helix_serial_number_;
}

PdbHelix::HelixResidueVector PdbHelix::GetHelixResidues()
{
    return helix_residues_;
}

PdbHelixClass PdbHelix::GetHelixClass()
{
    return helix_class_;
}

string PdbHelix::GetComment()
{
    return comment_;
}

double PdbHelix::GetHelixLength()
{
    return helix_length_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHelix::SetHelixId(const string helix_id)
{
    helix_id_ = helix_id;
}

void PdbHelix::SetHelixSerialNumber(int helix_serial_number)
{
    helix_serial_number_ = helix_serial_number;
}

void PdbHelix::SetHelixResidues(const HelixResidueVector helix_residues)
{
    helix_residues_.clear();
    for(PdbHelix::HelixResidueVector::const_iterator it = helix_residues.begin(); it != helix_residues.end(); it++)
    {
        helix_residues_.push_back(*it);
    }
}

void PdbHelix::AddHelixResidue(PdbHelixResidue *helix_residue)
{
    helix_residues_.push_back(helix_residue);
}

void PdbHelix::SetHelixClass(PdbHelixClass helix_class)
{
    helix_class_ = helix_class;
}

void PdbHelix::SetComment(const string &comment)
{
    comment_ = comment;
}

void PdbHelix::SetHelixLength(double helix_length)
{
    helix_length_ = helix_length;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

