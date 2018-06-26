#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbauthorsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbAuthorSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAuthorSection::PdbAuthorSection() : record_name_("AUTHOR"), author_(""){}

PdbAuthorSection::PdbAuthorSection(const std::string &record_name, const std::string &author)
{
    record_name_ = record_name;
    author_ = author;
}

PdbAuthorSection::PdbAuthorSection(std::stringstream& stream_block)
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
    this->SetAuthor( ss.str() );
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbAuthorSection::GetRecordName()
{
    return record_name_;
}

std::string PdbAuthorSection::GetAuthor()
{
    return author_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbAuthorSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbAuthorSection::SetAuthor(const std::string author)
{
    this->author_ = author;
	gmml::Trim( this->author_ );
	gmml::TrimSpaces( this->author_ );
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbAuthorSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << ", Author(s): " << author_ << std::endl << std::endl;
}
