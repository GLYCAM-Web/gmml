// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbheterogensynonym.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSynonym::PdbHeterogenSynonym() : heterogen_identifier_("") {}
PdbHeterogenSynonym::PdbHeterogenSynonym(const string &heterogen_identifier, const vector<string> &heterogen_synonyms) : heterogen_identifier_ (heterogen_identifier)
{
    heterogen_synonyms_.clear();
    for(vector<string>::const_iterator it = heterogen_synonyms.begin(); it != heterogen_synonyms.end(); it++)
    {
        heterogen_synonyms_.push_back(*it);
    }
}

PdbHeterogenSynonym::PdbHeterogenSynonym(stringstream& stream_block)
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
            is_heterogen_identifier_set = true;
        }

        ss << line.substr(15,55) << " ";

        getline(stream_block, line);
        temp = line;
    }
    heterogen_synonyms_ = Split(ss.str(), ";");
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenSynonym::GetHeterogenIdentifier()
{
    return heterogen_identifier_;
}

vector<string> PdbHeterogenSynonym::GetHeterogenSynonyms()
{
    return heterogen_synonyms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonym::SetHeterogenIdentifier(const string heterogen_identifier)
{
    heterogen_identifier_ = heterogen_identifier;
}

void PdbHeterogenSynonym::SetHeterogenSynonyms(const vector<string> heterogen_synonyms)
{
    heterogen_synonyms_.clear();
    for(vector<string>::const_iterator it = heterogen_synonyms.begin(); it != heterogen_synonyms.end(); it++)
    {
        heterogen_synonyms_.push_back(*it);
    }
}

void PdbHeterogenSynonym::AddHeterogenSynonym(const string heterogen_synonym)
{
    heterogen_synonyms_.push_back(heterogen_synonym);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonym::Print(ostream &out)
{
    out << "Heterogen ID: " << heterogen_identifier_ << "Heterogen Synonyms: ";
    for(vector<string>::iterator it = heterogen_synonyms_.begin(); it != heterogen_synonyms_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << endl;
}
