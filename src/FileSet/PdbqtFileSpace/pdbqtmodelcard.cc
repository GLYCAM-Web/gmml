#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbqtFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModelCard::PdbqtModelCard(){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

string PdbqtModelCard::GetRecordName(){
    return record_name_;
}

PdbqtModelCard::PdbqtModelMap PdbqtModelCard::GetModels(){
    return models_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbqtModelCard::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbqtModelCard::SetModels(PdbqtModelMap models){
    models_ = models;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtModelCard::Print(ostream &out)
{
}

