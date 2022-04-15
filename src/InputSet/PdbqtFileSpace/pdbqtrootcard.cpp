#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtRootCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtRootCard::PdbqtRootCard() : record_name_("ROOT")
{
    root_atoms_ = NULL;
}

PdbqtRootCard::PdbqtRootCard(std::ifstream &root_block, std::vector<PdbqtFileSpace::PdbqtAtomCard*>& ACV)
{
    root_atoms_ = NULL;
    std::string line;
    record_name_ = "ROOT";

    while (getline(root_block, line)){
        if(line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos)
        {
	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            root_block.seekg(offset, root_block.cur); //Go back one line
            root_atoms_ = new PdbqtFileSpace::PdbqtAtomCard(root_block);
	    ACV.push_back(root_atoms_);
            
        }

        else if(line.find("ENDROOT") != std::string::npos) //Quit Root section normally.
        {
	    break;
        }
	
	else { //The file stream has gone one line past the root section due to the absence of "ENDROOT". Should rewind by one line and quit RootCard.
	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            root_block.seekg(offset, root_block.cur); //Go back one line
	    break;
	}
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbqtRootCard::GetRecordName()
{
    return record_name_;
}
PdbqtFileSpace::PdbqtAtomCard* PdbqtRootCard::GetRootAtoms()
{
    return root_atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtRootCard::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbqtRootCard::SetRootAtoms(PdbqtFileSpace::PdbqtAtomCard* root_atoms)
{
    root_atoms_ = new PdbqtFileSpace::PdbqtAtomCard();
    root_atoms_->SetRecordName(root_atoms->GetRecordName());
    root_atoms_->SetAtoms(root_atoms->GetAtoms());
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtRootCard::Print(std::ostream &out)
{
    if(root_atoms_ != NULL)
        root_atoms_->Print(out);
}
