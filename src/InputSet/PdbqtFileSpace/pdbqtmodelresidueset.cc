#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtModelResidueSet;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModelResidueSet::PdbqtModelResidueSet()
{
    roots_ = NULL;
    atoms_ = NULL;
}

PdbqtModelResidueSet::PdbqtModelResidueSet(std::ifstream &residue_set_block)
{
    roots_ = NULL;
    atoms_ = NULL;
    branches_ = BranchCardVector();
    std::string line;

    //This vector contains atom cards from root and all downstream branches.They are then concatenated into a single atom card and stored in this object. 
    std::vector<PdbqtFileSpace::PdbqtAtomCard*> all_atom_cards; 

    while (getline(residue_set_block, line)){
        if(line.find("ROOT") != std::string::npos)
        {
            roots_ = new PdbqtFileSpace::PdbqtRootCard(residue_set_block, all_atom_cards);
        }

        else if(line.find("BRANCH") != std::string::npos)
        {
            std::vector<std::string> tokens = gmml::Split(line, " ");

	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            residue_set_block.seekg(offset, residue_set_block.cur); //Go back one line
            branches_.push_back(new PdbqtBranchCard(residue_set_block, all_atom_cards));
        }
	//If not ROOT,BRANCH,ENDROOT, or ENDBRANCH, the code has gone one line beyond the residue set section. So go back one line and exit. 
	else if (line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos){ 

	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            residue_set_block.seekg(offset, residue_set_block.cur); //Go back one line
	    //Get single atoms
	    all_atom_cards.push_back(new PdbqtFileSpace::PdbqtAtomCard(residue_set_block));
	}
        else if (line.find("TER") != std::string::npos){ //Do nothing with TER, ignore it
            continue;
        }
	else {
	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            residue_set_block.seekg(offset, residue_set_block.cur); //Go back one line
	    break;
	}
    }//while

    //TODO: Get all atoms
    atoms_ = new PdbqtFileSpace::PdbqtAtomCard();
    for (std::vector<PdbqtFileSpace::PdbqtAtomCard*>::iterator it = all_atom_cards.begin(); it != all_atom_cards.end(); it++){
	PdbqtAtomCard::PdbqtAtomMap card = (*it)->GetAtoms();
	for (PdbqtAtomCard::PdbqtAtomMap::iterator it2 = card.begin(); it2 != card.end(); it2++){
	    atoms_->AddAtom(it2->second);
	}
	
    }
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
PdbqtFileSpace::PdbqtRootCard* PdbqtModelResidueSet::GetRoots()
{
    return roots_;
}
PdbqtModelResidueSet::BranchCardVector PdbqtModelResidueSet::GetBranches()
{
    return branches_;
}
PdbqtFileSpace::PdbqtAtomCard* PdbqtModelResidueSet::GetAtoms()
{
    return atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtModelResidueSet::SetRoots(PdbqtFileSpace::PdbqtRootCard* roots)
{
    roots_ = roots;
}
void PdbqtModelResidueSet::SetBranches(BranchCardVector branches)
{
    branches_.clear();
    for(BranchCardVector::iterator it = branches.begin(); it != branches.end(); it++)
    {
        branches_.push_back(*it);
    }
}
void PdbqtModelResidueSet::AddBranch(PdbqtBranchCard* branch)
{
    branches_.push_back(branch);
}
void PdbqtModelResidueSet::SetAtoms(PdbqtFileSpace::PdbqtAtomCard *atoms)
{
    atoms_ = atoms;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtModelResidueSet::Print(std::ostream &out)
{
    out << "====================== Roots =====================" << std::endl;
    if(roots_ != NULL)
        roots_->Print(out);
    for(BranchCardVector::iterator it = branches_.begin(); it != branches_.end(); it++)
    {
        out << "====================== Main Branch =====================" << std::endl;
        (*it)->Print(out);
    }
}
