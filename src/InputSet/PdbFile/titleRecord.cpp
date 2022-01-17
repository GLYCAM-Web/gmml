#include <iostream>
#include <sstream>

#include "includes/InputSet/PdbFile/titleRecord.hpp"
#include "includes/utils.hpp"

using pdb::TitleRecord;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TitleRecord::TitleRecord(std::string record_name, std::string title)
{
    name_ = record_name;
    title_ = title;
}

TitleRecord::TitleRecord(std::stringstream& stream_block)
{
    std::string line;
    bool is_name_set = false;
    std::stringstream ss;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_name_set){
            name_ = line.substr(0,6);
            gmml::Trim(name_);
            is_name_set=true;
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
const std::string& TitleRecord::GetRecordName() const
{
    return name_;
}
const std::string& TitleRecord::GetTitle() const
{
    return title_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void TitleRecord::SetName(std::string record_name)
{
    name_ = record_name;
}

void TitleRecord::SetTitle(const std::string title)
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
void TitleRecord::Print(std::ostream &out) const
{
    out << "Record Name: " << name_ << ", Title: " << title_ << std::endl << std::endl;
}
