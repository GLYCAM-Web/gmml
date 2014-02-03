#include "../../../includes/FileSet/PdbFileSpace/pdbheadercard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

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

PdbHeaderCard::PdbHeaderCard(stringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        record_name_ = line.substr(0,6);
        Trim(record_name_);
        classification_ = line.substr(10,40);
        Trim(classification_);
        deposition_date_ = line.substr(50, 9);
        Trim(deposition_date_);
        identifier_code_ = line.substr(62,4);
        Trim(identifier_code_);
        getline(stream_block, line);
        temp = line;
    }
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
void PdbHeaderCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << ", Classification: " << classification_ <<
           ", Deposition Date: " << deposition_date_ << ", Identifier Code: " << identifier_code_ << endl << endl;
}
