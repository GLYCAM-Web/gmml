#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbcaveatsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbCaveatSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCaveatSection::PdbCaveatSection() : record_name_("CAVEAT"), caveat_(""){}

PdbCaveatSection::PdbCaveatSection(const std::string &record_name, const std::string &caveat)
{
    record_name_ = record_name;
    caveat_ = caveat;
}

PdbCaveatSection::PdbCaveatSection(std::stringstream& stream_block)
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
        ss << line.substr(20,70);

        getline(stream_block, line);
        temp = line;
    }
    caveat_ = ss.str();
    // gmml::Trim(caveat_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string PdbCaveatSection::GetRecordName()
{
    return record_name_;
}

std::string PdbCaveatSection::GetCaveat()
{
    return caveat_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbCaveatSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbCaveatSection::SetCaveat(const std::string caveat)
{
    caveat_ = caveat;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbCaveatSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << ", Caveat: " << caveat_ << std::endl << std::endl;
}
