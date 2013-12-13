#include "../../../includes/FileSet/PdbFileSpace/pdbheadercard.h"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeaderCard::PdbHeaderCard() : record_name_("HEADER"), classification_(""), deposition_date_(""), identifier_code_(""){}

PdbHeaderCard::PdbHeaderCard(const string &record_name, const string &classification, const string &deposition_date, const string &identifier_code)
{
    record_name_ = record_name;
    classification_ = classification;
    deposition_date_ = deposition_date;
    identifier_code_ = identifier_code;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbHeaderCard::GetRecordName()
{
    return record_name_;
}

string PdbHeaderCard::GetClassification()
{
    return classification_;
}

string PdbHeaderCard::GetDepositionDate()
{
    return deposition_date_;
}

string PdbHeaderCard::GetIdentifierCode()
{
    return identifier_code_;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
void PdbHeaderCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbHeaderCard::SetClassification(const string classification)
{
    classification_ = classification;
}

void PdbHeaderCard::SetDepositionDate(const string deposition_date)
{
    deposition_date_ = deposition_date;
}

void PdbHeaderCard::SetIdentificationCode(const string identifier_code){
    identifier_code_ = identifier_code;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
