#include "../../../includes/FileSet/PdbFileSpace/pdbnummodelcard.hpp"

using namespace std;
using namespace PdbFileSpace;


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbNumModelCard::PdbNumModelCard() : record_name_("NUMMDL"), number_of_models_(0) {}
PdbNumModelCard::PdbNumModelCard(const string &record_name, int number_of_models) : record_name_(record_name), number_of_models_(number_of_models) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbNumModelCard::GetRecordName()
{
    return record_name_;
}

int PdbNumModelCard::GetNumberOfModels()
{
    return number_of_models_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbNumModelCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbNumModelCard::SetNumberOfModels(int number_of_models)
{
    number_of_models_ = number_of_models;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
