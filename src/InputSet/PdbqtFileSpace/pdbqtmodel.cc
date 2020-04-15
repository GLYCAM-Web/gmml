#include "../../../includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"

using PdbqtFileSpace::PdbqtModel;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModel::PdbqtModel() {}

PdbqtModel::PdbqtModel(std::ifstream &model_block)
{
    std::string line;
    std::stringstream residue_set_block;
    remarks_ = RemarkCardVector();
    torsional_dof_cards_ = TorsionalDoFCardVector();
    model_compound_card_ = NULL;

    //This is the MODEL line if file contains MODEL at all. If not, then the file doesn't contain MODEL, could be ROOT,BRANCH etc.
    //Should the latter be the case, go back one line since the current line contains useful information. It needs to be read again below.   
    getline(model_block, line); 
    if (line.find("MODEL") != std::string::npos){ //If there is MODEL entry
            std::vector<std::string> tokens = gmml::Split(line, " ");
            if(tokens.size() == 1 || tokens.at(1).empty())
                model_serial_number_ = gmml::iNotSet;
            else
                model_serial_number_ = gmml::ConvertString<int>(tokens.at(1));
        }
    else { //If there isn't a MODEL entry, the entire file is one model. Also go up one line
        model_serial_number_ = 1;
	int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
        model_block.seekg(offset, model_block.cur);//Go back one line
    }


    while(getline(model_block,line)){
	//MODEL
        if (line.find("MODEL") != std::string::npos){ //If there is MODEL entry
	    //MODEL has been dealt with above. This will be duplicate MODEL lines, ignore. 
        }

        // COMPND
        else if(line.find("COMPND") != std::string::npos)
        {
            model_compound_card_ = new PdbqtFileSpace::PdbqtCompoundCard(line);
            //getline(model_block, line);
        }
        // REMARK
        else if(line.find("REMARK") != std::string::npos)
        {
            PdbqtRemarkCard* remark = new PdbqtRemarkCard(line);
            remarks_.push_back(remark);
            //getline(model_block, line);
        }
        // ROOT/ENDROOT/ATOM/HETATM/BRANCH/ENDBRANCH
        else if(line.find("ROOT") != std::string::npos || line.find("ATOM") != std::string::npos
              || line.find("HETATM") != std::string::npos || line.find("ENDROOT") != std::string::npos
              || line.find("BRANCH") != std::string::npos || line.find("ENDBRANCH") != std::string::npos)
        {
	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            model_block.seekg(offset, model_block.cur);//Go back one line
            model_residue_set_ = new PdbqtFileSpace::PdbqtModelResidueSet(model_block);
        }
        // TORSDOF
        else if (line.find("TORSDOF") != std::string::npos)
        {
            torsional_dof_cards_.push_back(new PdbqtTorsionalDoFCard(line));
            //getline(model_block, line);
        }
	else if (line.find("ENDMDL") != std::string::npos){ //Exit upon ENDMDL is the normal outcome.
	    break;
	}
        else if (line.find("TER") != std::string::npos){
            continue;
        }
	else{  //If not the keywords above, the file stream has gone one line beyond the model block. Go back one line and exit
	    //int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            //model_block.seekg(offset, model_block.cur);//Go back one line
	    //break;
            throw PdbqtFileProcessingException(__LINE__, "Illegal record detected in model section");
        }
    }//while
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int PdbqtModel::GetModelSerialNumber()
{
    return model_serial_number_;
}
PdbqtFileSpace::PdbqtCompoundCard* PdbqtModel::GetModelCompoundCard()
{
    return model_compound_card_;
}
PdbqtFileSpace::PdbqtModelResidueSet* PdbqtModel::GetModelResidueSet()
{
    return model_residue_set_;
}
PdbqtModel::RemarkCardVector PdbqtModel::GetRemarks()
{
    return remarks_;
}
PdbqtModel::TorsionalDoFCardVector PdbqtModel::GetTorsionalDoFCards()
{
    return torsional_dof_cards_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbqtModel::SetModelSerialNumber(int model_serial_number)
{
    model_serial_number_ = model_serial_number;
}
void PdbqtModel::SetModelCompundCard(PdbqtFileSpace::PdbqtCompoundCard* model_compound_card)
{
    model_compound_card_ = model_compound_card;
}
void PdbqtModel::SetModelResidueSet(PdbqtFileSpace::PdbqtModelResidueSet* model_residue_set)
{
    model_residue_set_ = model_residue_set;
}
void PdbqtModel::SetRemarks(RemarkCardVector remarks)
{
    remarks_.clear();
    for(RemarkCardVector::iterator it = remarks.begin(); it != remarks.end(); it++)
    {
        remarks_.push_back(*it);
    }
}
void PdbqtModel::AddRemark(PdbqtRemarkCard *remark)
{
    remarks_.push_back(remark);
}
void PdbqtModel::SetTorsionalDofCards(TorsionalDoFCardVector torsional_dof_cards)
{
    torsional_dof_cards_.clear();
    for(TorsionalDoFCardVector::iterator it = torsional_dof_cards.begin(); it != torsional_dof_cards.end(); it++)
    {
        torsional_dof_cards_.push_back(*it);
    }
}
void PdbqtModel::AddTorsionalDoFCard(PdbqtTorsionalDoFCard *torsional_dof_cards)
{
    torsional_dof_cards_.push_back(torsional_dof_cards);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtModel::Print(std::ostream &out)
{
    out << "====================== Compound Card =====================" << std::endl;
    model_compound_card_->Print(out);
    out << "====================== Remarks =====================" << std::endl;
    for(RemarkCardVector::iterator it = remarks_.begin(); it != remarks_.end(); it++)
        (*it)->Print(out);
    out << std::endl
        << "====================== Residue Set =====================" << std::endl;
    model_residue_set_->Print(out);
    out << std::endl
        << "====================== Torsional DOF =====================" << std::endl;
    out << std::endl;
    for(TorsionalDoFCardVector::iterator it = torsional_dof_cards_.begin(); it != torsional_dof_cards_.end(); it++)
        (*it)->Print(out);
}
