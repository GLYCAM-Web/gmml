// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbHeterogenSynonymCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard() : heterogen_identifier_("") {}
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard(const std::string &heterogen_identifier, const std::vector<std::string> &heterogen_synonyms) : heterogen_identifier_ (heterogen_identifier)
{
    heterogen_synonyms_.clear();
    for(std::vector<std::string>::const_iterator it = heterogen_synonyms.begin(); it != heterogen_synonyms.end(); it++)
    {
        heterogen_synonyms_.push_back(*it);
    }
}

PdbHeterogenSynonymCard::PdbHeterogenSynonymCard(std::stringstream& stream_block)
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
    std::string temp_synonym = ss.str();
    heterogen_synonyms_ = gmml::Split(gmml::Trim(temp_synonym), ";");
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHeterogenSynonymCard::GetHeterogenIdentifier()
{
    return heterogen_identifier_;
}

std::vector<std::string> PdbHeterogenSynonymCard::GetHeterogenSynonymCards()
{
    return heterogen_synonyms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymCard::SetHeterogenIdentifier(const std::string heterogen_identifier)
{
    heterogen_identifier_ = heterogen_identifier;
}

void PdbHeterogenSynonymCard::SetHeterogenSynonyms(const std::vector<std::string> heterogen_synonyms)
{
    heterogen_synonyms_.clear();
    for(std::vector<std::string>::const_iterator it = heterogen_synonyms.begin(); it != heterogen_synonyms.end(); it++)
    {
        heterogen_synonyms_.push_back(*it);
    }
}

void PdbHeterogenSynonymCard::AddHeterogenSynonym(const std::string heterogen_synonym)
{
    heterogen_synonyms_.push_back(heterogen_synonym);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymCard::Print(std::ostream &out)
{
    out << "Heterogen ID: " << heterogen_identifier_ << "Heterogen Synonyms: ";
    for(std::vector<std::string>::iterator it = heterogen_synonyms_.begin(); it != heterogen_synonyms_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << std::endl;
}
