#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbcaveatsection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCaveatSection::PdbCaveatSection() : record_name_("CAVEAT"), caveat_(""){}

PdbCaveatSection::PdbCaveatSection(const string &record_name, const string &caveat)
{
    record_name_ = record_name;
    caveat_ = caveat;
}

PdbCaveatSection::PdbCaveatSection(stringstream& stream_block)
{
    string line;
    bool is_record_name_set = false;
    stringstream ss;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
        ss << line.substr(20,70);

        getline(stream_block, line);
        temp = line;
    }
    caveat_ = ss.str();
    // Trim(caveat_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbCaveatSection::GetRecordName()
{
    return record_name_;
}

string PdbCaveatSection::GetCaveat()
{
    return caveat_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbCaveatSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbCaveatSection::SetCaveat(const string caveat)
{
    caveat_ = caveat;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbCaveatSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << ", Caveat: " << caveat_ << endl << endl;
}
