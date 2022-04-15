#include "../../../includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbNumModelCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbNumModelCard::PdbNumModelCard() : record_name_("NUMMDL"), number_of_models_(0) {}
PdbNumModelCard::PdbNumModelCard(const std::string &record_name, int number_of_models) : record_name_(record_name), number_of_models_(number_of_models) {}

PdbNumModelCard::PdbNumModelCard(std::stringstream& stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        record_name_ = line.substr(0,6);
        gmml::Trim(record_name_);
        if(line.substr(10,4) == "    ")
            number_of_models_ = gmml::iNotSet;
        else
            number_of_models_ = gmml::ConvertString<int>(line.substr(10,4));

        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbNumModelCard::GetRecordName()
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
void PdbNumModelCard::SetRecordName(const std::string record_name)
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
void PdbNumModelCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_
        << ", Number of Models: ";
    if(number_of_models_ != gmml::iNotSet)
        out << number_of_models_;
    else
        out << " ";
    out << std::endl << std::endl;
}
