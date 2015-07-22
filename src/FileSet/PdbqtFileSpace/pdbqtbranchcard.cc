#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtBranchCard::PdbqtBranchCard() : record_name_("BRANCH"){}

PdbqtBranchCard::PdbqtBranchCard(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    if(line.find("BRANCH") != string::npos)
    {
        vector<string> tokens = Split(line, " ");
        record_name_ = tokens.at(0);
        solid_atom_serial_number_ = ConvertString<int>(tokens.at(1));
        rotatable_atom_serial_number_ = ConvertString<int>(tokens.at(2));
        sub_branches_ = BranchCardVector();
        getline(stream_block, line);
        stringstream rotatable_set_block;
        while(line.find("ATOM") != string::npos || line.find("HETATM") != string::npos)
        {
            rotatable_set_block << line << endl;
            getline(stream_block, line);
        }
        rotatable_atom_set_ = new PdbqtAtomCard(rotatable_set_block);
        stringstream end_branch_line;
        end_branch_line << "ENDBRANCH " << setw(3) << solid_atom_serial_number_ << " " << setw(3) << rotatable_atom_serial_number_;
        while(line.find(end_branch_line.str()) == string::npos)
        {
            if(line.find("BRANCH") != string::npos)
            {
                vector<string> subbranch_tokens = Split(line, " ");
                int solid_atom_serial_number_of_subbranch = ConvertString<int>(subbranch_tokens.at(1));
                int rotatable_atom_serial_number_of_subbranch = ConvertString<int>(subbranch_tokens.at(2));
                stringstream end_subbranch_line;
                end_subbranch_line << "ENDBRANCH " << setw(3) << solid_atom_serial_number_of_subbranch << " "
                                   << setw(3) << rotatable_atom_serial_number_of_subbranch;
                stringstream subbranch_block;
                subbranch_block << line << endl;
                getline(stream_block, line);
                while(line.find(end_subbranch_line.str()) == string::npos)
                {
                    subbranch_block << line << endl;
                    getline(stream_block, line);
                }
                if(line.find(end_subbranch_line.str()) != string::npos)
                {
                    subbranch_block << line << endl;
                    sub_branches_.push_back(new PdbqtBranchCard(subbranch_block));
                    getline(stream_block, line);
                }
            }
        }
    }
}

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
PdbqtAtomCard* PdbqtBranchCard::GetRotatableAtomSet()
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
void PdbqtBranchCard::SetRotatableAtomSet(PdbqtAtomCard* rotatable_atom_set)
{
    rotatable_atom_set_ = new PdbqtAtomCard();
    rotatable_atom_set_->SetRecordName(rotatable_atom_set->GetRecordName());
    rotatable_atom_set_->SetAtoms(rotatable_atom_set->GetAtoms());
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
    out << "BRANCH " << solid_atom_serial_number_ << " " << rotatable_atom_serial_number_ << endl;
    out << "=========================== Rotatable =======================" << endl;
    rotatable_atom_set_->Print(out);
    out << endl << "============================== Sub Branches =============================" << endl;
    for(BranchCardVector::iterator it = sub_branches_.begin(); it != sub_branches_.end(); it++)
        (*it)->Print(out);
}

