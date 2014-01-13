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
    stringstream ss;
    bool is_record_name_set = false;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }

        ss << line;
        PdbHelix* helix = new PdbHelix(ss);
        helixes_[line.substr(11,3)] = helix;
        getline(stream_block, line);
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

