#include "../../../includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbNumModelCard::PdbNumModelCard() : record_name_("NUMMDL"), number_of_models_(0) {}
PdbNumModelCard::PdbNumModelCard(const string &record_name, int number_of_models) : record_name_(record_name), number_of_models_(number_of_models) {}

PdbNumModelCard::PdbNumModelCard(stringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        record_name_ = line.substr(0,6);
        Trim(record_name_);
        if(line.substr(10,4) == "    ")
            number_of_models_ = iNotSet;
        else
            number_of_models_ = ConvertString<int>(line.substr(10,4));

        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbNumModelCard::GetRecordName()
{
    return record_name_;
}

int PdbNumModelCard::GetNumberOfModels()
{
    return number_of_models_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbNumModelCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbNumModelCard::SetNumberOfModels(int number_of_models)
{
    number_of_models_ = number_of_models;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbNumModelCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_
        << ", Number of Models: ";
    if(number_of_models_ != iNotSet)
        out << number_of_models_;
    else
        out << " ";
    out << endl << endl;
}
