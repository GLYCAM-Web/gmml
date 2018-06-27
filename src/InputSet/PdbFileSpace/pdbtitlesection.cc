#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbTitleSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbTitleSection::PdbTitleSection() : record_name_("TITLE"), title_(""){}

PdbTitleSection::PdbTitleSection(const std::string &record_name, const std::string &title)
{
    record_name_ = record_name;
    title_ = title;
}

PdbTitleSection::PdbTitleSection(std::stringstream& stream_block)
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
    this->SetTitle( ss.str() );
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbTitleSection::GetRecordName()
{
    return record_name_;
}

std::string PdbTitleSection::GetTitle()
{
    return title_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbTitleSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbTitleSection::SetTitle(const std::string title)
{
    this->title_ = title;
	gmml::Trim( this->title_ );
	gmml::TrimSpaces( this->title_ );
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbTitleSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << ", Title: " << title_ << std::endl << std::endl;
}
