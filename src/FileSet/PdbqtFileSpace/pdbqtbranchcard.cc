#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtBranchCard::PdbqtBranchCard(){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbqtBranchCard::GetRecordName()
{
    return record_name_;
}
int PdbqtBranchCard::GetSolidAtomSerialNumber()
{
    return solid_atom_serial_number_;
}
int PdbqtBranchCard::GetRotatbleAtomSerialNumber()
{
    return rotatable_atom_serial_number_;
}
PdbqtBranchCard::AtomCardVector PdbqtBranchCard::GetRotatableAtomSet()
{
    return rotatable_atom_set_;
}
PdbqtBranchCard::BranchCardVector PdbqtBranchCard::GetSubBranches()
{
    return sub_branches_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtBranchCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}
void PdbqtBranchCard::SetSolidAtomSerialNumber(int solid_atom_serial_number)
{
    solid_atom_serial_number_ = solid_atom_serial_number;
}
void PdbqtBranchCard::SetRotatableAtomSerialNumber(int rotatable_atom_serial_number)
{
    rotatable_atom_serial_number_ = rotatable_atom_serial_number;
}
void PdbqtBranchCard::SetRotatableAtomSet(AtomCardVector rotatable_atom_set)
{
    rotatable_atom_set_.clear();
    for(AtomCardVector::iterator it = rotatable_atom_set.begin(); it != rotatable_atom_set.end(); it++)
    {
        rotatable_atom_set_.push_back(*it);
    }
}
void PdbqtBranchCard::AddRotatableAtom(PdbqtAtomCard* rotatable_atom)
{
    rotatable_atom_set_.push_back(rotatable_atom);
}
void PdbqtBranchCard::SetSubBranches(BranchCardVector sub_branch)
{
    sub_branches_.clear();
    for(BranchCardVector::iterator it = sub_branch.begin(); it != sub_branch.end(); it++)
    {
        sub_branches_.push_back(*it);
    }
}
void PdbqtBranchCard::AddSubBranch(PdbqtBranchCard* sub_branch)
{
    sub_branches_.push_back(sub_branch);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtBranchCard::Print(ostream &out)
{
}

