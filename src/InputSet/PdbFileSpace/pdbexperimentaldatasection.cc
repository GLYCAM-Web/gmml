#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbexperimentaldatasection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbExperimentalDataSection::PdbExperimentalDataSection() : record_name_("EXPDTA"), experimental_data_(""){}

PdbExperimentalDataSection::PdbExperimentalDataSection(const string &record_name, const string &experimental_data)
{
    record_name_ = record_name;
    experimental_data_ = experimental_data;
}

PdbExperimentalDataSection::PdbExperimentalDataSection(stringstream& stream_block)
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
        ss << line.substr(10,70);

        getline(stream_block, line);
        temp = line;
    }
    experimental_data_ = ss.str();
    Trim(experimental_data_);
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbExperimentalDataSection::GetRecordName()
{
    return record_name_;
}

string PdbExperimentalDataSection::GetExperimentalData()
{
    return experimental_data_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbExperimentalDataSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbExperimentalDataSection::SetExperimentalData(const string experimental_data)
{
    experimental_data_ = experimental_data;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbExperimentalDataSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << ", Experimental Data: " << experimental_data_ << endl << endl;
}
