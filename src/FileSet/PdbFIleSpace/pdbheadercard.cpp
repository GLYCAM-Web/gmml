#include "../../../includes/FileSet/PdbFileSpace/pdbheadercard.h"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeaderClass::PdbHeaderClass() : record_name_("HEADER"), classification_(""), deposition_date_(""), identifier_code_(""){}

PdbHeaderClass::PdbHeaderClass(string record_name, string classification, string deposition_date, string identifier_code)
{
    record_name_ = record_name;
    classification_ = classification;
    deposition_date_ = deposition_date;
    identifier_code_ = identifier_code;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbHeaderClass::GetRecordName()
{
    return record_name_;
}

string PdbHeaderClass::GetClassification()
{
    return classification_;
}

string PdbHeaderClass::GetDepositionDate()
{
    return deposition_date_;
}

string PdbHeaderClass::GetIdentifierCode()
{
    return identifier_code_;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
void PdbHeaderClass::SetRecordName(string record_name)
{
    record_name_ = record_name;
}

void PdbHeaderClass::SetClassification(string classification)
{
    classification_ = classification;
}

void PdbHeaderClass::SetDepositionDate(string deposition_date)
{
    deposition_date_ = deposition_date;
}

void PdbHeaderClass::SetIdentificationCode(string identifier_code){
    identifier_code_ = identifier_code;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
