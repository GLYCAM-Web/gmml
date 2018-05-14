#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbexperimentaldatasection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbExperimentalDataSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbExperimentalDataSection::PdbExperimentalDataSection() : record_name_("EXPDTA"), experimental_data_(""){}

PdbExperimentalDataSection::PdbExperimentalDataSection(const std::string &record_name, const std::string &experimental_data)
{
    record_name_ = record_name;
    experimental_data_ = experimental_data;
}

PdbExperimentalDataSection::PdbExperimentalDataSection(std::stringstream& stream_block)
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
    experimental_data_ = ss.str();
    gmml::Trim(experimental_data_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbExperimentalDataSection::GetRecordName()
{
    return record_name_;
}

std::string PdbExperimentalDataSection::GetExperimentalData()
{
    return experimental_data_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbExperimentalDataSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbExperimentalDataSection::SetExperimentalData(const std::string experimental_data)
{
    experimental_data_ = experimental_data;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbExperimentalDataSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << ", Experimental Data: " << experimental_data_ << std::endl << std::endl;
}
