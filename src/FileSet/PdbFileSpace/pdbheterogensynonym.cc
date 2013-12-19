// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbheterogensynonym.hpp"

using namespace std;
using namespace PdbFileSpace;

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

