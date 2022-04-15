#include "../../../includes/InputSet/PdbFileSpace/pdbsourcecard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbSourceCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSourceCard::PdbSourceCard() {}
PdbSourceCard::PdbSourceCard(std::string &line)
{

    // gmml::Trim(record_name_);
    std::string temp = line;
    if (!gmml::Trim(temp).empty())
    {
      record_name_ = line.substr(0, 6);
      std::size_t position = line.find(":");
      if (position!=std::string::npos)
         {
           token_ = line.substr(10,position-10);
           // gmml::Trim(token_);
           value_ = line.substr(position+1, 78-position);
           // gmml::Trim(value_);
         }
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string PdbSourceCard::GetRecordName(){
    return record_name_;
}

std::string PdbSourceCard::GetToken(){
    return token_;
}

std::string PdbSourceCard::GetValue(){
    return value_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSourceCard::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

void PdbSourceCard::SetToken(const std::string token){
    token_ = token;
}

void PdbSourceCard::SetValue(const std::string value){
    value_ = value;
}


//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbSourceCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "Token:" << token_;
    out << "Value:" << value_;
    out << std::endl;
}
