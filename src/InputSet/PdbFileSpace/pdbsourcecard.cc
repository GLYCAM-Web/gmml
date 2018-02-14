#include "../../../includes/InputSet/PdbFileSpace/pdbsourcecard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSourceCard::PdbSourceCard() {}
PdbSourceCard::PdbSourceCard(string &line)
{
    record_name_ = line.substr(0, 6);
    Trim(record_name_);

    if (!Trim(line).empty())
    {
      std::size_t position = line.find(":");
      if (position!=std::string::npos)
         {
           token_ = line.substr(10,position-10);
           Trim(token_);
           value_ = line.substr(position+1, 78-position);
           Trim(value_);
         }
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbSourceCard::GetRecordName(){
    return record_name_;
}

string PdbSourceCard::GetToken(){
    return token_;
}

string PdbSourceCard::GetValue(){
    return value_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSourceCard::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbSourceCard::SetToken(const string token){
    token_ = token;
}

void PdbSourceCard::SetValue(const string value){
    value_ = value;
}


//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbSourceCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "Token:" << token_;
    out << "Value:" << value_;
    out << endl;
}
