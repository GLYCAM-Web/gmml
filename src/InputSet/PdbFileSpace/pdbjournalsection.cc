#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbjournalsection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbJournalSection::PdbJournalSection() /*: record_name_("JRNL"), authors_(""),
                                         title_(""), editors_(""), reference_(""),
                                         publisher_(""), reference_nums_(""),
                                         pmid_(""), doi_("")*/{}

PdbJournalSection::PdbJournalSection(stringstream& stream_block)
{
    string line;
    bool is_record_name_set = false;
    bool is_title_started = false;
    bool is_reference_started = false;
    bool is_publisher_started = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set)
        {
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
        string subrecord = line.substr(12,4);
        Trim(subrecord);
        if(subrecord == "AUTH")
        {
          std::size_t start_position = 19;
          std::size_t end_position = line.find(",");
          while (end_position!=std::string::npos)
             {
               string new_author = line.substr(start_position,end_position);
               authors_.push_back(new_author);
               start_position = end_position;
               end_position = line.find(",",start_position+1);
             }
        }
        else if(subrecord == "TITL")
        {
          if(!is_title_started)
          {
            title_ = line.substr(19,59);
            Trim(title_);
            is_title_started = true;
          }
          else
          {
            string titl = line.substr(19,59);
            Trim(titl);
            title_.append(" ");
            title_.append(titl);
            Trim(title_);
          }
        }
        else if(subrecord == "EDIT")
        {
          std::size_t start_position = 19;
          std::size_t end_position = line.find(",");
          while (end_position!=std::string::npos)
             {
               string new_editor = line.substr(start_position,end_position);
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
            Trim(reference_);
            is_reference_started = true;
          }
          else
          {
            string ref = line.substr(19,59);
            Trim(ref);
            reference_.append(" ");
            reference_.append(ref);
            Trim(reference_);
          }
        }
        else if(subrecord == "PUBL")
        {
          if(!is_publisher_started)
          {
            publisher_ = line.substr(19,59);
            Trim(publisher_);
            is_publisher_started = true;
          }
          else
          {
            string publ = line.substr(19,59);
            Trim(publ);
            publisher_.append(" ");
            publisher_.append(publ);
          }
        }
        else if(subrecord == "REFN")
        {
          string refn = line.substr(35,29);
          Trim(refn);
          reference_nums_.push_back(refn);
        }
        else if(subrecord == "PMID")
        {
            string new_pmid = line.substr(19,59);
            Trim(new_pmid);
            pmid_.append(new_pmid);
        }
        else if(subrecord == "DOI")
        {
          string new_doi = line.substr(19,59);
          Trim(new_doi);
          doi_.append(new_doi);
        }
    }
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbJournalSection::GetRecordName()
{
    return record_name_;
}

std::vector<string> PdbJournalSection::GetAuthors()
{
    return authors_;
}

string PdbJournalSection::GetTitle()
{
    return title_;
}

std::vector<string> PdbJournalSection::GetEditors()
{
    return editors_;
}

string PdbJournalSection::GetReference()
{
    return reference_;
}

string PdbJournalSection::GetPublisher()
{
    return publisher_;
}

std::vector<string> PdbJournalSection::GetReferenceNumbers()
{
    return reference_nums_;
}

string PdbJournalSection::GetPMID()
{
    return pmid_;
}

string PdbJournalSection::GetDOI()
{
    return doi_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbJournalSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbJournalSection::SetAuthors(std::vector<string> authors)
{
    authors_ = authors;
}

void PdbJournalSection::SetTitle(const string title)
{
    title_ = title;
}

void PdbJournalSection::SetEditors(std::vector<string> editors)
{
    editors_ = editors;
}

void PdbJournalSection::SetReference(const string reference)
{
    reference_ = reference;
}

void PdbJournalSection::SetPublisher(const string publisher)
{
    publisher_ = publisher;
}

void PdbJournalSection::SetReferenceNumbers(std::vector<string> reference_nums)
{
    reference_nums_ = reference_nums;
}

void PdbJournalSection::SetPMID(const string pmid)
{
    pmid_ = pmid;
}

void PdbJournalSection::SetDOI(const string doi)
{
    doi_ = doi;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbJournalSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl;
    out << ", Authors: ";
    for(vector<string>::iterator it = authors_.begin(); it != authors_.end(); it++)
    {
      out << *it << " ";
    }
    out << endl;
    out << "Title: " << title_ << endl;
    out << "Editors: ";
    for(vector<string>::iterator it = editors_.begin(); it != editors_.end(); it++)
    {
      out << *it << " ";
    }
    out << endl;
    out << "Reference: " << reference_ << endl;
    out << "Publisher: " << publisher_ << endl;
    out << "Reference Numbers: ";
    for(vector<string>::iterator it = reference_nums_.begin();
        it != reference_nums_.end(); it++)
    {
      out << *it << " ";
    }
    out << endl;
    out << "PMID: " << pmid_ << endl;
    out << "DOI: " << doi_ << endl;
}
