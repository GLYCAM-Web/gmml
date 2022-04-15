// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbHeterogenNameCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenNameCard::PdbHeterogenNameCard() : heterogen_identifier_(""), heterogen_name_("") {}
PdbHeterogenNameCard::PdbHeterogenNameCard(const std::string &heterogen_identifier, const std::string &heterogen_name) :
    heterogen_identifier_(heterogen_identifier), heterogen_name_(heterogen_name) {}

PdbHeterogenNameCard::PdbHeterogenNameCard(std::stringstream& stream_block)
{
    std::string line;
    bool is_heterogen_identifier_set = false;
    std::stringstream ss;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_heterogen_identifier_set){
            heterogen_identifier_ = line.substr(11,3);
            gmml::Trim(heterogen_identifier_);
            is_heterogen_identifier_set = true;
        }

        ss << line.substr(15,55) << " ";

        getline(stream_block, line);
        temp = line;
    }
    heterogen_name_ = ss.str();
    heterogen_name_ = gmml::Trim(heterogen_name_);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHeterogenNameCard::GetHeterogenIdentifier()
{
    return heterogen_identifier_;
}

std::string PdbHeterogenNameCard::GetHeterogenName()
{
    return heterogen_name_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenNameCard::SetHeterogenIdentifier(const std::string heterogen_identifier)
{
    heterogen_identifier_ = heterogen_identifier;
}

void PdbHeterogenNameCard::SetHeterogenName(const std::string heterogen_name)
{
    heterogen_name_ = heterogen_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenNameCard::Print(std::ostream &out)
{
    out << "Heterogen ID: " << heterogen_identifier_ << ", Heterogen Name: " << heterogen_name_ << std::endl;
}
