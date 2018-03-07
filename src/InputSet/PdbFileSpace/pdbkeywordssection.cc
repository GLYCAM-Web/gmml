#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbkeywordssection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbKeywordsSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbKeywordsSection::PdbKeywordsSection() : record_name_("KEYWDS"), keywords_(""){}

PdbKeywordsSection::PdbKeywordsSection(const std::string &record_name, const std::string &keywords)
{
    record_name_ = record_name;
    keywords_ = keywords;
}

PdbKeywordsSection::PdbKeywordsSection(std::stringstream& stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    std::stringstream ss;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
        ss << line.substr(10,70);

        getline(stream_block, line);
        temp = line;
    }
    keywords_ = ss.str();
    gmml::Trim(keywords_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbKeywordsSection::GetRecordName()
{
    return record_name_;
}

std::string PdbKeywordsSection::GetKeywords()
{
    return keywords_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbKeywordsSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbKeywordsSection::SetKeywords(const std::string keywords)
{
    keywords_ = keywords;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbKeywordsSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << ", Keywords: " << keywords_ << std::endl << std::endl;
}
