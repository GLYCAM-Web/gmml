// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard() : heterogen_identifier_("") {}
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard(const string &heterogen_identifier, const vector<string> &heterogen_synonyms) : heterogen_identifier_ (heterogen_identifier)
{
    heterogen_synonyms_.clear();
    for(vector<string>::const_iterator it = heterogen_synonyms.begin(); it != heterogen_synonyms.end(); it++)
    {
        heterogen_synonyms_.push_back(*it);
    }
}

PdbHeterogenSynonymCard::PdbHeterogenSynonymCard(stringstream& stream_block)
{
    string line;
    bool is_heterogen_identifier_set = false;
    stringstream ss;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_heterogen_identifier_set){
            heterogen_identifier_ = line.substr(11,3);
            Trim(heterogen_identifier_);
            is_heterogen_identifier_set = true;
        }

        ss << line.substr(15,55) << " ";

        getline(stream_block, line);
        temp = line;
    }
    string temp_synonym = ss.str();
    heterogen_synonyms_ = Split(Trim(temp_synonym), ";");
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenSynonymCard::GetHeterogenIdentifier()
{
    return heterogen_identifier_;
}

vector<string> PdbHeterogenSynonymCard::GetHeterogenSynonymCards()
{
    return heterogen_synonyms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymCard::SetHeterogenIdentifier(const string heterogen_identifier)
{
    heterogen_identifier_ = heterogen_identifier;
}

void PdbHeterogenSynonymCard::SetHeterogenSynonyms(const vector<string> heterogen_synonyms)
{
    heterogen_synonyms_.clear();
    for(vector<string>::const_iterator it = heterogen_synonyms.begin(); it != heterogen_synonyms.end(); it++)
    {
        heterogen_synonyms_.push_back(*it);
    }
}

void PdbHeterogenSynonymCard::AddHeterogenSynonym(const string heterogen_synonym)
{
    heterogen_synonyms_.push_back(heterogen_synonym);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymCard::Print(ostream &out)
{
    out << "Heterogen ID: " << heterogen_identifier_ << "Heterogen Synonyms: ";
    for(vector<string>::iterator it = heterogen_synonyms_.begin(); it != heterogen_synonyms_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << endl;
}
