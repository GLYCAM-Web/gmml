#include <iostream>

#include "includes/InputSet/PdbFile/authorRecord.hpp"
#include "includes/utils.hpp"

using pdb::AuthorRecord;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
AuthorRecord::AuthorRecord() : record_name_("AUTHOR"), author_(""){}

AuthorRecord::AuthorRecord(std::string record_name, std::string author)
{
    record_name_ = record_name;
    author_ = author;
}

AuthorRecord::AuthorRecord(std::stringstream& stream_block)
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
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void AuthorRecord::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void AuthorRecord::SetAuthor(const std::string author)
{
    this->author_ = author;
	gmml::Trim( this->author_ );
	gmml::TrimSpaces( this->author_ );
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void AuthorRecord::Print(std::ostream &out) const
{
    out << "Record Name: " << this->GetRecordName() << ", Author(s): " << this->GetAuthor() << std::endl << std::endl;
}
