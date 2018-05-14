// Author Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbhelixcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbhelixsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbHelixSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHelixSection::PdbHelixSection() : record_name_("HELIX") {}

PdbHelixSection::PdbHelixSection(const std::string &record_name) : record_name_(record_name) {}

PdbHelixSection::PdbHelixSection(std::stringstream& stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }

        std::stringstream ss;
        ss << line << std::endl;
        PdbHelixCard* helix = new PdbHelixCard(ss);
        helix_cards_[helix->GetHelixId()] = helix;
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHelixSection::GetRecordName()
{
    return record_name_;
}

PdbHelixSection::HelixCardMap PdbHelixSection::GetHelixCards()
{
    return helix_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHelixSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHelixSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "============= Helixes ===========" << std::endl;
    for(PdbHelixSection::HelixCardMap::iterator it = helix_cards_.begin(); it != helix_cards_.end(); it++)
    {
        out << "Helix ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
