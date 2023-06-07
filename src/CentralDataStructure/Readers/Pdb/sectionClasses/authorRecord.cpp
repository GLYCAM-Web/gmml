#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/authorRecord.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <iostream>
#include <iomanip> //setw
#include <cmath> //ceil

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
    while (!codeUtils::Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            codeUtils::Trim(record_name_);
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
    codeUtils::Trim( this->author_ );
    codeUtils::removeMultipleSpaces( this->author_ );
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void AuthorRecord::Print(std::ostream &out) const
{
    out << "Record Name: " << this->GetRecordName() << ", Author(s): " << this->GetAuthor() << std::endl << std::endl;
}

void AuthorRecord::Write(std::ostream& stream) const
{
    const int MAX_AUTHOR_LENGTH_IN_LINE = 70;
    stream << std::left << std::setw(6) << this->GetRecordName()
                     << std::left << std::setw(2) << " ";
    if((int)this->GetAuthor().length() > MAX_AUTHOR_LENGTH_IN_LINE)
    {
        stream << std::left << std::setw(70) << this->GetAuthor().substr(0,MAX_AUTHOR_LENGTH_IN_LINE)
                         << std::endl;

        int counter = ceil((double)(this->GetAuthor().length()) / MAX_AUTHOR_LENGTH_IN_LINE);
        for(int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << std::left << std::setw(6) << this->GetRecordName()
                                 << std::left << std::setw(2) << " "
                                 << std::right << std::setw(2) << i
                                 << std::left << std::setw(70) << this->GetAuthor().substr(MAX_AUTHOR_LENGTH_IN_LINE*(i-1), MAX_AUTHOR_LENGTH_IN_LINE)
                                 << std::endl;
            }
            else
            {
                stream << std::left << std::setw(6) << this->GetRecordName()
                                 << std::left << std::setw(2) << " "
                                 << std::right << std::setw(2) << i
                                 << std::left << std::setw(70) << this->GetAuthor().substr(MAX_AUTHOR_LENGTH_IN_LINE*(i-1), this->GetAuthor().length()-MAX_AUTHOR_LENGTH_IN_LINE*(i-1))
                                 << std::endl;
            }
        }
    }
    else
    {
        stream << std::right << std::setw(2) << " "
                << std::left << std::setw(70) << this->GetAuthor()
                << std::endl;
    }
}

