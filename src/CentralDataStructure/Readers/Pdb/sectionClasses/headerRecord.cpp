#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/headerRecord.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <iomanip> //setw

using pdb::HeaderRecord;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
HeaderRecord::HeaderRecord() : record_name_("HEADER"), classification_(" "), deposition_date_(" "), identifier_code_(" "){}

HeaderRecord::HeaderRecord(const std::string &record_name, const std::string &classification, const std::string &deposition_date, const std::string &identifier_code)
{
    record_name_ = record_name;
    classification_ = classification;
    deposition_date_ = deposition_date;
    identifier_code_ = identifier_code;
}

HeaderRecord::HeaderRecord(std::stringstream& stream_block)
{
  record_name_ = "";
  classification_ = "";
  deposition_date_ = "";
  identifier_code_ = "";
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!codeUtils::Trim(temp).empty())
    {
        //record_name_ = line.substr(0,6);
        //gmml::Trim(record_name_);
        //classification_ = line.substr(10,40);
        //gmml::Trim(classification_);
        //deposition_date_ = line.substr(50, 9);
        //gmml::Trim(deposition_date_);
        //identifier_code_ = line.substr(62,4);
        //gmml::Trim(identifier_code_);
        this->SetRecordName(codeUtils::RemoveWhiteSpace(line.substr(0, 6)));
        this->SetClassification(codeUtils::RemoveWhiteSpace(line.substr(10, 40)));
        this->SetDepositionDate(codeUtils::RemoveWhiteSpace(line.substr(50, 9)));
        this->SetIdentificationCode(codeUtils::RemoveWhiteSpace(line.substr(62, 4)));
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string HeaderRecord::GetRecordName() const
{
    return record_name_;
}

std::string HeaderRecord::GetClassification() const
{
    return classification_;
}

std::string HeaderRecord::GetDepositionDate() const
{
    return deposition_date_;
}

std::string HeaderRecord::GetIdentifierCode() const
{
    return identifier_code_;
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
void HeaderRecord::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void HeaderRecord::SetClassification(const std::string classification)
{
    classification_ = classification;
}

void HeaderRecord::SetDepositionDate(const std::string deposition_date)
{
    deposition_date_ = deposition_date;
}

void HeaderRecord::SetIdentificationCode(const std::string identifier_code)
{
    identifier_code_ = identifier_code;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void HeaderRecord::Print(std::ostream &out) const
{
    out << "Record Name: " << this->GetRecordName() << ", Classification: " << this->GetClassification() <<
           ", Deposition Date: " << this->GetDepositionDate() << ", Identifier Code: " << this->GetIdentifierCode() << std::endl << std::endl;
}


void HeaderRecord::Write(std::ostream& stream) const
{
    stream << std::left << std::setw(6) << this->GetRecordName()
           << std::left << std::setw(4) << " "
           << std::left << std::setw(40) << this->GetClassification()
           << std::left << std::setw(9) << this->GetDepositionDate()
           << std::left << std::setw(3) << " "
           << std::right << std::setw(4) << this->GetIdentifierCode()
           << std::left << std::setw(14) << " "
           << std::endl;
}
