// Author Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbhelixcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbhelixsection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHelixSection::PdbHelixSection() : record_name_("HELIX") {}

PdbHelixSection::PdbHelixSection(const string &record_name) : record_name_(record_name) {}

PdbHelixSection::PdbHelixSection(stringstream& stream_block)
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

        stringstream ss;
        ss << line << endl;
        PdbHelixCard* helix = new PdbHelixCard(ss);
        helix_cards_[helix->GetHelixId()] = helix;
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHelixSection::GetRecordName()
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
void PdbHelixSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHelixSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "============= Helixes ===========" << endl;
    for(PdbHelixSection::HelixCardMap::iterator it = helix_cards_.begin(); it != helix_cards_.end(); it++)
    {
        out << "Helix ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
