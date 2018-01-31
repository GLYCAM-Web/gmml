// Created by: Dave Montgomery

#include "../../../includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"


using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbRemarkSection::PdbRemarkSection() : record_name_("REMARK"){}

PdbRemarkSection::PdbRemarkSection(const std::string& record_name, const std::string& remark_cards)
{
    record_name_ = record_name;
    remark_cards_ = remark_cards;
}

PdbRemarkSection::PdbRemarkSection(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(line.find("REMARK") != string::npos)
        {
            if(!is_record_name_set)
            {
                record_name_ = line.substr(0,6);
                Trim(record_name_);
                is_record_name_set=true;
            }
            stringstream remark_cardstream;
            while(line.find("REMARK") != string::npos)
            {
                remark_cardstream << line << endl;
                getline(stream_block,line);
                temp = line;
            }
            remark_cards_ = remark_cardstream.str();
        }

    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

string PdbRemarkSection::GetRecordName(){
    return record_name_;
}

string PdbRemarkSection::GetRemarks(){
    return remark_cards_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbRemarkSection::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbRemarkSection::SetRemarks(const string remark_cards){
    remark_cards_ = remark_cards;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbRemarkSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           remark_cards_ << endl;

    // for(PdbRemarkSection::PdbRemarkMap::iterator it = remark_cards_.begin(); it != remark_cards_.end(); it++)
    // {
    //     out << "Remark Serial Number: ";
    //     if((it)->first != iNotSet)
    //         out << (it)->first << endl;
    //     else
    //         out << " " << endl;
    //     (it)->second->Print();
    //     out << endl;
    // }
}
