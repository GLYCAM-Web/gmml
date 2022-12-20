#include "includes/CentralDataStructure/Readers/Pdb/SectionClasses/journalRecord.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <iomanip> //setw
#include <iostream>
#include <sstream>
#include <cmath> // ceil

using pdb::JournalRecord;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
JournalRecord::JournalRecord(){};
JournalRecord::JournalRecord(std::stringstream& stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    bool is_title_started = false;
    bool is_reference_started = false;
    bool is_publisher_started = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!codeUtils::Trim(temp).empty())
    {
        if(!is_record_name_set)
        {
            record_name_ = line.substr(0,6);
            codeUtils::Trim(record_name_);
            is_record_name_set=true;
        }
        text_.append(line.substr(12, 67));
        std::string subrecord = line.substr(12,4);
        codeUtils::Trim(subrecord);
        if(subrecord == "AUTH")
        {
          std::size_t start_position = 19;
          std::size_t end_position = line.find(",");
          if (end_position!=std::string::npos)
             {
               std::string new_author = line.substr(start_position,end_position-start_position);
               authors_.push_back(new_author);
               start_position = end_position;
               end_position = line.find(",",start_position+1);
             }
          else
          {
            std::string new_author = line.substr(start_position,79-start_position);
            authors_.push_back(new_author);
          }
        }
        else if(subrecord == "TITL")
        {
          if(!is_title_started)
          {
            title_ = line.substr(19,59);
            codeUtils::Trim(title_);
            is_title_started = true;
          }
          else
          {
            std::string titl = line.substr(19,59);
            codeUtils::Trim(titl);
            title_.append(" ");
            title_.append(titl);
            codeUtils::Trim(title_);
          }
        }
        else if(subrecord == "EDIT")
        {
          std::size_t start_position = 19;
          std::size_t end_position = line.find(",");
          while (end_position!=std::string::npos)
             {
               std::string new_editor = line.substr(start_position,end_position-start_position);
               editors_.push_back(new_editor);
               start_position = end_position;
               end_position = line.find(",",start_position+1);
             }
        }
        else if(subrecord == "REF")
        {
          if(!is_reference_started)
          {
            reference_ = line.substr(19,59);
            codeUtils::Trim(reference_);
            is_reference_started = true;
          }
          else
          {
            std::string ref = line.substr(19,59);
            codeUtils::Trim(ref);
            reference_.append(" ");
            reference_.append(ref);
            codeUtils::Trim(reference_);
          }
        }
        else if(subrecord == "PUBL")
        {
          if(!is_publisher_started)
          {
            publisher_ = line.substr(19,59);
            codeUtils::Trim(publisher_);
            is_publisher_started = true;
          }
          else
          {
            std::string publ = line.substr(19,59);
            codeUtils::Trim(publ);
            publisher_.append(" ");
            publisher_.append(publ);
          }
        }
        else if(subrecord == "REFN")
        {
          std::string refn = line.substr(35,29);
          codeUtils::Trim(refn);
          reference_nums_.push_back(refn);
        }
        else if(subrecord == "PMID")
        {
            std::string new_pmid = line.substr(19,59);
            codeUtils::Trim(new_pmid);
            pmid_.append(new_pmid);
        }
        else if(subrecord == "DOI")
        {
          std::string new_doi = line.substr(19,59);
          this->SetDOI( new_doi );
        }
    getline(stream_block, line);
    temp = line;
    }
}
//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void JournalRecord::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void JournalRecord::SetAuthors(std::vector<std::string> authors)
{
    authors_ = authors;
}

void JournalRecord::SetTitle(const std::string title)
{
    title_ = title;
}

void JournalRecord::SetEditors(std::vector<std::string> editors)
{
    editors_ = editors;
}

void JournalRecord::SetReference(const std::string reference)
{
    reference_ = reference;
}

void JournalRecord::SetPublisher(const std::string publisher)
{
    publisher_ = publisher;
}

void JournalRecord::SetReferenceNumbers(std::vector<std::string> reference_nums)
{
    reference_nums_ = reference_nums;
}

void JournalRecord::SetPMID(const std::string pmid)
{
    pmid_ = pmid;
}

void JournalRecord::SetDOI(const std::string doi)
{
    this->doi_ += doi;
	codeUtils::Trim( this->doi_ );
    codeUtils::removeMultipleSpaces( this->doi_ );
	// I know this is weird to Trim spaces then add one space back, but it was/is necessary. DT
	this->doi_ += " ";
}

void JournalRecord::SetText(const std::string text)
{
    text_ = text;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void JournalRecord::Print(std::ostream &out) const
{
    out << "Record Name: " << record_name_ << std::endl;
    out << "Authors: ";
    for(auto &author : this->GetAuthors())
    {
      out << author << " ";
    }
    out << std::endl;
    out << "Title: " << title_ << std::endl;
    out << "Editors: ";
    for(auto &editor : this->GetEditors())
    {
        out << editor << " ";
    }
    out << std::endl;
    out << "Reference: " << reference_ << std::endl;
    out << "Publisher: " << publisher_ << std::endl;
    out << "Reference Numbers: ";
    for(auto &refNum : this->GetReferenceNumbers())
    {
        out << refNum << " ";
    }
    out << std::endl;
    out << "PMID: " << pmid_ << std::endl;
    out << "DOI: " << doi_ << std::endl;
}

void JournalRecord::Write(std::ostream& stream) const
{
    const int MAX_JOURNAL_LENGTH_IN_LINE = 67;
    stream << std::left << std::setw(6) << this->GetRecordName()
                     << std::left << std::setw(6) << " ";
    if((int)this->GetText().length() > MAX_JOURNAL_LENGTH_IN_LINE)
    {
        stream << std::left << std::setw(67) << this->GetText().substr(0,MAX_JOURNAL_LENGTH_IN_LINE)
                         << std::endl;

        int counter = ceil((double)(this->GetText().length()) / MAX_JOURNAL_LENGTH_IN_LINE);
        for(int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << std::left << std::setw(6) << this->GetRecordName()
                                 << std::left << std::setw(6) << " "
                                 << std::left << std::setw(67) << this->GetText().substr(MAX_JOURNAL_LENGTH_IN_LINE*(i-1), MAX_JOURNAL_LENGTH_IN_LINE)
                                 << std::endl;
            }
            else
            {
                stream << std::left << std::setw(6) << this->GetRecordName()
                                 << std::left << std::setw(6) << " "
                                 << std::left << std::setw(67) << this->GetText().substr(MAX_JOURNAL_LENGTH_IN_LINE*(i-1), this->GetText().length()-MAX_JOURNAL_LENGTH_IN_LINE*(i-1))
                                 << std::endl;
            }
        }
    }
    else
    {
        stream << std::right << std::setw(6) << " "
                << std::left << std::setw(67) << this
                ->GetText()
                << std::endl;
    }
}

