#include "../../../includes/InputSet/PdbFileSpace/pdbmodeltypesection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbModelTypeSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModelTypeSection::PdbModelTypeSection() : record_name_("MDLTYP") {}

PdbModelTypeSection::PdbModelTypeSection(const std::string &record_name, const std::vector<std::string> &comments) : record_name_(record_name), comments_(comments) {}

PdbModelTypeSection::PdbModelTypeSection(std::stringstream& stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    std::stringstream ss;
    std::string temp_comments;
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
    temp_comments = ss.str();
    comments_ = gmml::Split(temp_comments, ",");
    for(std::vector<std::string>::iterator it = comments_.begin(); it != comments_.end(); it++)
    {
        gmml::Trim(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbModelTypeSection::GetRecordName()
{
    return record_name_;
}

std::vector<std::string> PdbModelTypeSection::GetComments()
{
    return comments_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbModelTypeSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbModelTypeSection::SetComments(const std::vector<std::string> comments)
{
    comments_.clear();
    for(std::vector<std::string>::const_iterator it = comments.begin(); it != comments.end(); it++)
    {
        comments_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbModelTypeSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl << "Comments: ";
    for(std::vector<std::string>::iterator it = comments_.begin(); it != comments_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << std::endl << std::endl;
}
