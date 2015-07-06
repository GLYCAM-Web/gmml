#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
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

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int PdbqtModel::GetModelSerialNumber()
{
    return model_serial_number_;
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
}

