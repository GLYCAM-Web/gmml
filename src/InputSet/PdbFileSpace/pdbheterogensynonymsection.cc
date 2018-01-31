// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymsection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSynonymSection::PdbHeterogenSynonymSection() : record_name_("HETSYN") {}
PdbHeterogenSynonymSection::PdbHeterogenSynonymSection(const string& record_name) : record_name_(record_name) {}
PdbHeterogenSynonymSection::PdbHeterogenSynonymSection(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
        stringstream heterogen_synonym_block;
        heterogen_synonym_block << line << endl;
        string heterogen_identifier = line.substr(11,3);

        getline(stream_block, line);
        temp = line;

        while (!Trim(temp).empty() && line.substr(11,3) == heterogen_identifier){
            heterogen_synonym_block << line << endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbHeterogenSynonymCard* heterogen_synonym = new PdbHeterogenSynonymCard(heterogen_synonym_block);
        heterogen_identifier = Trim(heterogen_identifier);
        heterogens_synonym_cards_[heterogen_identifier] = heterogen_synonym;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenSynonymSection::GetRecordName()
{
    return record_name_;
}

PdbHeterogenSynonymSection::HeterogenSynonymCardMap PdbHeterogenSynonymSection::GetHeterogensSynonymCards()
{
    return heterogens_synonym_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "============ Heterogen Synonyms ===========" << endl;
    for(PdbHeterogenSynonymSection::HeterogenSynonymCardMap::iterator it = heterogens_synonym_cards_.begin(); it != heterogens_synonym_cards_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
