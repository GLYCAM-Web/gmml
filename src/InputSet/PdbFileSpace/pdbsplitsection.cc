#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbsplitsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbSplitSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSplitSection::PdbSplitSection() : record_name_("SPLIT"), split_(""){}

PdbSplitSection::PdbSplitSection(const std::string &record_name, const std::string &split)
{
    record_name_ = record_name;
    split_ = split;
}

PdbSplitSection::PdbSplitSection(std::stringstream& stream_block)
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
    split_ = ss.str();
    gmml::Trim(split_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbSplitSection::GetRecordName()
{
    return record_name_;
}

std::string PdbSplitSection::GetSplit()
{
    return split_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbSplitSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbSplitSection::SetSplit(const std::string split)
{
    split_ = split;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSplitSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << ", Split: " << split_ << std::endl << std::endl;
}
