#include "../../../includes/FileSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbqtFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModel::PdbqtModel() {}

PdbqtModel::PdbqtModel(stringstream &model_block)
{
    string line;
    stringstream residue_set_block;
    remarks_ = RemarkCardVector();
    torsional_dof_cards_ = TorsionalDoFCardVector();
    model_compound_card_ = NULL;
    getline(model_block, line);
    if(line.find("MODEL") != string::npos)
    {
        vector<string> tokens = Split(line, " ");
        if(tokens.size() == 1 || tokens.at(1).empty())
            model_serial_number_ = iNotSet;
        else
            model_serial_number_ = ConvertString<int>(tokens.at(1));
        getline(model_block,line);
        // COMPND
        if(line.find("COMPND") != string::npos)
        {
            model_compound_card_ = new PdbqtCompoundCard(line);
            getline(model_block, line);
        }
        // REMARK
        while(line.find("REMARK") != string::npos)
        {
            PdbqtRemarkCard* remark = new PdbqtRemarkCard(line);
            remarks_.push_back(remark);
            getline(model_block, line);
        }
        // ROOT/ENDROOT/ATOM/HETATM/BRANCH/ENDBRANCH
        while(line.find("ROOT") != string::npos || line.find("ATOM") != string::npos
              || line.find("HETATM") != string::npos || line.find("ENDROOT") != string::npos
              || line.find("BRANCH") != string::npos || line.find("ENDBRANCH") != string::npos)
        {
            residue_set_block << line << endl;
            getline(model_block, line);
        }
        model_residue_set_ = new PdbqtModelResidueSet(residue_set_block);
        // TORSDOF
        while(line.find("TORSDOF") != string::npos)
        {
            torsional_dof_cards_.push_back(new PdbqtTorsionalDoFCard(line));
            getline(model_block, line);
        }
    }
    else
    {
        model_serial_number_ = 1;
        string temp = line;
        // COMPND
        if(line.find("COMPND") != string::npos)
        {
            model_compound_card_ = new PdbqtCompoundCard(line);
            getline(model_block, line);
        }
        // REMARK
        while(line.find("REMARK") != string::npos)
        {
            PdbqtRemarkCard* remark = new PdbqtRemarkCard(line);
            remarks_.push_back(remark);
            getline(model_block, line);
        }
        // ROOT/ENDROOT/ATOM/HETATM/BRANCH/ENDBRANCH
        while(line.find("ATOM") != string::npos || line.find("ANISOU") != string::npos
              || line.find("TER") != string::npos || line.find("HETATM") != string::npos)
        {
            residue_set_block << line << endl;
            getline(model_block, line);
        }
        model_residue_set_ = new PdbqtModelResidueSet(residue_set_block);
        // TORSDOF
        while(line.find("TORSDOF") != string::npos)
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
PdbqtCompoundCard* PdbqtModel::GetModelCompoundCard()
{
    return model_compound_card_;
}
PdbqtModelResidueSet* PdbqtModel::GetModelResidueSet()
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
void PdbqtModel::SetModelCompundCard(PdbqtCompoundCard* model_compound_card)
{
    model_compound_card_ = model_compound_card;
}
void PdbqtModel::SetModelResidueSet(PdbqtModelResidueSet* model_residue_set)
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
void PdbqtModel::Print(ostream &out)
{
    out << "====================== Compound Card =====================" << endl;
    model_compound_card_->Print(out);
    out << "====================== Remarks =====================" << endl;
    for(RemarkCardVector::iterator it = remarks_.begin(); it != remarks_.end(); it++)
        (*it)->Print(out);
    out << endl
        << "====================== Residue Set =====================" << endl;
    model_residue_set_->Print(out);
    out << endl
        << "====================== Torsional DOF =====================" << endl;
    out << endl;
    for(TorsionalDoFCardVector::iterator it = torsional_dof_cards_.begin(); it != torsional_dof_cards_.end(); it++)
        (*it)->Print(out);
}

