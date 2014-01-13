// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogensynonym.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard() : record_name_("HETSYN") {}
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard(const string& record_name) : record_name_(record_name) {}
PdbHeterogenSynonymCard::PdbHeterogenSynonymCard(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }
        stringstream heterogen_synonym_block;
        heterogen_synonym_block << line << endl;
        char heterogen_identifier = line.substr(11,3);

        getline(stream_block, line);

        while (!Trim(line).empty() && line.substr(11,3) == heterogen_identifier){
            heterogen_synonym_block << line << endl;
            getline(stream_block, line);
        }
        PdbHeterogenSynonym* heterogen_synonym = new PdbHeterogenSynonym(heterogen_synonym_block);
        heterogens_synonyms_[heterogen_identifier] = heterogen_synonym;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenSynonymCard::GetRecordName()
{
    return record_name_;
}

PdbHeterogenSynonymCard::HeterogenSynonymMap PdbHeterogenSynonymCard::GetHeterogensSynonyms()
{
    return heterogens_synonyms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenSynonymCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////


