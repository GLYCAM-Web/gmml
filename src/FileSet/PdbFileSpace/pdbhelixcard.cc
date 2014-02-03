// Author Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbhelixcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbhelix.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHelixCard::PdbHelixCard() : record_name_("HELIX") {}

PdbHelixCard::PdbHelixCard(const string &record_name) : record_name_(record_name) {}

PdbHelixCard::PdbHelixCard(stringstream& stream_block)
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
        PdbHelix* helix = new PdbHelix(ss);
        helixes_[helix->GetHelixId()] = helix;
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHelixCard::GetRecordName()
{
    return record_name_;
}

PdbHelixCard::HelixMap PdbHelixCard::GetHelixes()
{
    return helixes_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHelixCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHelixCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "============= Helixes ===========" << endl;
    for(PdbHelixCard::HelixMap::iterator it = helixes_.begin(); it != helixes_.end(); it++)
    {
        out << "Helix ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
