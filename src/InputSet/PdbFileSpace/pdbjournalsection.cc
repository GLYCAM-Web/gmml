#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbjournalsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbJournalSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbJournalSection::PdbJournalSection() /*: record_name_("JRNL"), authors_(""),
                                         title_(""), editors_(""), reference_(""),
                                         publisher_(""), reference_nums_(""),
                                         pmid_(""), doi_("")*/{}

PdbJournalSection::PdbJournalSection(std::stringstream& stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    bool is_title_started = false;
    bool is_reference_started = false;
    bool is_publisher_started = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set)
        {
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
        text_.append(line.substr(12, 67));
        std::string subrecord = line.substr(12,4);
        gmml::Trim(subrecord);
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
            gmml::Trim(title_);
            is_title_started = true;
          }
          else
          {
            std::string titl = line.substr(19,59);
            gmml::Trim(titl);
            title_.append(" ");
            title_.append(titl);
            gmml::Trim(title_);
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
            gmml::Trim(reference_);
            is_reference_started = true;
          }
          else
          {
            std::string ref = line.substr(19,59);
            gmml::Trim(ref);
            reference_.append(" ");
            reference_.append(ref);
            gmml::Trim(reference_);
          }
        }
        else if(subrecord == "PUBL")
        {
          if(!is_publisher_started)
          {
            publisher_ = line.substr(19,59);
            gmml::Trim(publisher_);
            is_publisher_started = true;
          }
          else
          {
            std::string publ = line.substr(19,59);
            gmml::Trim(publ);
            publisher_.append(" ");
            publisher_.append(publ);
          }
        }
        else if(subrecord == "REFN")
        {
          std::string refn = line.substr(35,29);
          gmml::Trim(refn);
          reference_nums_.push_back(refn);
        }
        else if(subrecord == "PMID")
        {
            std::string new_pmid = line.substr(19,59);
            gmml::Trim(new_pmid);
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
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbJournalSection::GetRecordName()
{
    return record_name_;
}

std::vector<std::string> PdbJournalSection::GetAuthors()
{
    return authors_;
}

std::string PdbJournalSection::GetTitle()
{
    return title_;
}

std::vector<std::string> PdbJournalSection::GetEditors()
{
    return editors_;
}

std::string PdbJournalSection::GetReference()
{
    return reference_;
}

std::string PdbJournalSection::GetPublisher()
{
    return publisher_;
}

std::vector<std::string> PdbJournalSection::GetReferenceNumbers()
{
    return reference_nums_;
}

std::string PdbJournalSection::GetPMID()
{
    return pmid_;
}

std::string PdbJournalSection::GetDOI()
{
    return doi_;
}

std::string PdbJournalSection::GetText()
{
    return text_;
}
//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbJournalSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbJournalSection::SetAuthors(std::vector<std::string> authors)
{
    authors_ = authors;
}

void PdbJournalSection::SetTitle(const std::string title)
{
    title_ = title;
}

void PdbJournalSection::SetEditors(std::vector<std::string> editors)
{
    editors_ = editors;
}

void PdbJournalSection::SetReference(const std::string reference)
{
    reference_ = reference;
}

void PdbJournalSection::SetPublisher(const std::string publisher)
{
    publisher_ = publisher;
}

void PdbJournalSection::SetReferenceNumbers(std::vector<std::string> reference_nums)
{
    reference_nums_ = reference_nums;
}

void PdbJournalSection::SetPMID(const std::string pmid)
{
    pmid_ = pmid;
}

void PdbJournalSection::SetDOI(const std::string doi)
{
    this->doi_ += doi;
	gmml::Trim( this->doi_ );
    gmml::TrimSpaces( this->doi_ );
	// I know this is weird to Trim spaces then add one space back, but it was/is necessary. DT
	this->doi_ += " ";
}

void PdbJournalSection::SetText(const std::string text)
{
    text_ = text;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbJournalSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl;
    out << "Authors: ";
    for(std::vector<std::string>::iterator it = authors_.begin(); it != authors_.end(); it++)
    {
      out << *it << " ";
    }
    out << std::endl;
    out << "Title: " << title_ << std::endl;
    out << "Editors: ";
    for(std::vector<std::string>::iterator it = editors_.begin(); it != editors_.end(); it++)
    {
      out << *it << " ";
    }
    out << std::endl;
    out << "Reference: " << reference_ << std::endl;
    out << "Publisher: " << publisher_ << std::endl;
    out << "Reference Numbers: ";
    for(std::vector<std::string>::iterator it = reference_nums_.begin();
        it != reference_nums_.end(); it++)
    {
      out << *it << " ";
    }
    out << std::endl;
    out << "PMID: " << pmid_ << std::endl;
    out << "DOI: " << doi_ << std::endl;
}
