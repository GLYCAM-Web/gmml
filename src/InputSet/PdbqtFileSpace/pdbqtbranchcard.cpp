#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtBranchCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtBranchCard::PdbqtBranchCard() : record_name_("BRANCH"){}

PdbqtBranchCard::PdbqtBranchCard(std::ifstream &stream_block, std::vector<PdbqtFileSpace::PdbqtAtomCard*>& ACV)
{
    sub_branches_ = BranchCardVector();
    std::string line;
    bool subbranch_exists = false;
    while (getline(stream_block, line)){
        if(line.find("BRANCH") != std::string::npos)
        {
	    if (!subbranch_exists){
                std::vector<std::string> tokens = gmml::Split(line, " ");
                record_name_ = tokens.at(0);
                solid_atom_serial_number_ = gmml::ConvertString<int>(tokens.at(1));
                rotatable_atom_serial_number_ = gmml::ConvertString<int>(tokens.at(2));

	        while (getline(stream_block, line)){
                    if(line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos)
                    {
		        int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            	        stream_block.seekg(offset, stream_block.cur); //Go back one line
                        rotatable_atom_set_ = new PdbqtFileSpace::PdbqtAtomCard(stream_block);
			ACV.push_back(rotatable_atom_set_);
                    }
		    else{ //Either BRANCH or ENDBRANCH
		        int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            	        stream_block.seekg(offset, stream_block.cur); //Go back one line
		        if (line.find("ENDBRANCH") != std::string::npos){
			    subbranch_exists = false;
		        }
		        else if (line.find("BRANCH") != std::string::npos){
			    subbranch_exists = true;
		        }
		        break;
		    }
	        }
	    }//if subbranch doesn't exist. 

	    else{//If there is a subbranch
	        int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            	stream_block.seekg(offset, stream_block.cur); //Go back one line
                sub_branches_.push_back(new PdbqtBranchCard(stream_block, ACV));
	    }
        }
	else if (line.find("ENDBRANCH") != std::string::npos){
	    break;
	}
	else {//If current line is not BRANCH,ATOM,HETATM or ENDBRANCH, the code has gone one line beyond the Branch section. Go back one line and exit.
	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            stream_block.seekg(offset, stream_block.cur); //Go back one line
	    break;
	}
    }//while
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbqtBranchCard::GetRecordName()
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
PdbqtFileSpace::PdbqtAtomCard* PdbqtBranchCard::GetRotatableAtomSet()
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
void PdbqtBranchCard::SetRecordName(const std::string record_name)
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
void PdbqtBranchCard::SetRotatableAtomSet(PdbqtFileSpace::PdbqtAtomCard* rotatable_atom_set)
{
    rotatable_atom_set_ = new PdbqtFileSpace::PdbqtAtomCard();
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
void PdbqtBranchCard::Print(std::ostream &out)
{
    out << "BRANCH " << solid_atom_serial_number_ << " " << rotatable_atom_serial_number_ << std::endl;
    out << "=========================== Rotatable =======================" << std::endl;
    rotatable_atom_set_->Print(out);
    out << std::endl << "============================== Sub Branches =============================" << std::endl;
    for(BranchCardVector::iterator it = sub_branches_.begin(); it != sub_branches_.end(); it++)
        (*it)->Print(out);
}
