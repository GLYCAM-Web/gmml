#include "../../../includes/InputSet/PdbFileSpace/pdbheadercard.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbHeaderCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeaderCard::PdbHeaderCard() : record_name_("HEADER"), classification_(" "), deposition_date_(" "), identifier_code_(" "){}

PdbHeaderCard::PdbHeaderCard(const std::string &record_name, const std::string &classification, const std::string &deposition_date, const std::string &identifier_code)
{
    record_name_ = record_name;
    classification_ = classification;
    deposition_date_ = deposition_date;
    identifier_code_ = identifier_code;
}

PdbHeaderCard::PdbHeaderCard(std::stringstream& stream_block)
{
  record_name_ = "";
  classification_ = "";
  deposition_date_ = "";
  identifier_code_ = "";
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        record_name_ = line.substr(0,6);
        gmml::Trim(record_name_);
        classification_ = line.substr(10,40);
        gmml::Trim(classification_);
        deposition_date_ = line.substr(50, 9);
        gmml::Trim(deposition_date_);
        identifier_code_ = line.substr(62,4);
        gmml::Trim(identifier_code_);
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbHeaderCard::GetRecordName()
{
    return record_name_;
}

std::string PdbHeaderCard::GetClassification()
{
    return classification_;
}

std::string PdbHeaderCard::GetDepositionDate()
{
    return deposition_date_;
}

std::string PdbHeaderCard::GetIdentifierCode()
{
    return identifier_code_;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
void PdbHeaderCard::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbHeaderCard::SetClassification(const std::string classification)
{
    classification_ = classification;
}

void PdbHeaderCard::SetDepositionDate(const std::string deposition_date)
{
    deposition_date_ = deposition_date;
}

void PdbHeaderCard::SetIdentificationCode(const std::string identifier_code){
    identifier_code_ = identifier_code;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeaderCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << ", Classification: " << classification_ <<
           ", Deposition Date: " << deposition_date_ << ", Identifier Code: " << identifier_code_ << std::endl << std::endl;
}
