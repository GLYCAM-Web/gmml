#include "../../../includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtModel;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModel::PdbqtModel() {}

PdbqtModel::PdbqtModel(std::stringstream &model_block)
{
    std::string line;
    std::stringstream residue_set_block;
    remarks_ = RemarkCardVector();
    torsional_dof_cards_ = TorsionalDoFCardVector();
    model_compound_card_ = NULL;
    getline(model_block, line);
    if(line.find("MODEL") != std::string::npos)
    {
        std::vector<std::string> tokens = gmml::Split(line, " ");
        if(tokens.size() == 1 || tokens.at(1).empty())
            model_serial_number_ = gmml::iNotSet;
        else
            model_serial_number_ = gmml::ConvertString<int>(tokens.at(1));
        getline(model_block,line);
        // COMPND
        if(line.find("COMPND") != std::string::npos)
        {
            model_compound_card_ = new PdbqtFileSpace::PdbqtCompoundCard(line);
            getline(model_block, line);
        }
        // REMARK
        while(line.find("REMARK") != std::string::npos)
        {
            PdbqtRemarkCard* remark = new PdbqtRemarkCard(line);
            remarks_.push_back(remark);
            getline(model_block, line);
        }
        // ROOT/ENDROOT/ATOM/HETATM/BRANCH/ENDBRANCH
        while(line.find("ROOT") != std::string::npos || line.find("ATOM") != std::string::npos
              || line.find("HETATM") != std::string::npos || line.find("ENDROOT") != std::string::npos
              || line.find("BRANCH") != std::string::npos || line.find("ENDBRANCH") != std::string::npos)
        {
            residue_set_block << line << std::endl;
            getline(model_block, line);
        }
        model_residue_set_ = new PdbqtFileSpace::PdbqtModelResidueSet(residue_set_block);
        // TORSDOF
        while(line.find("TORSDOF") != std::string::npos)
        {
            torsional_dof_cards_.push_back(new PdbqtTorsionalDoFCard(line));
            getline(model_block, line);
        }
    }
    else
    {
        model_serial_number_ = 1;
        std::string temp = line;
        // COMPND
        if(line.find("COMPND") != std::string::npos)
        {
            model_compound_card_ = new PdbqtFileSpace::PdbqtCompoundCard(line);
            getline(model_block, line);
        }
        // REMARK
        while(line.find("REMARK") != std::string::npos)
        {
            PdbqtRemarkCard* remark = new PdbqtRemarkCard(line);
            remarks_.push_back(remark);
            getline(model_block, line);
        }
        // ROOT/ENDROOT/ATOM/HETATM/BRANCH/ENDBRANCH
        while(line.find("ATOM") != std::string::npos || line.find("ANISOU") != std::string::npos
              || line.find("TER") != std::string::npos || line.find("HETATM") != std::string::npos)
        {
            residue_set_block << line << std::endl;
            getline(model_block, line);
        }
        model_residue_set_ = new PdbqtFileSpace::PdbqtModelResidueSet(residue_set_block);
        // TORSDOF
        while(line.find("TORSDOF") != std::string::npos)
        {
            torsional_dof_cards_.push_back(new PdbqtTorsionalDoFCard(line));
            getline(model_block, line);
        }
    }
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
