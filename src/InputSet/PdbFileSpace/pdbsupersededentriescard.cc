#include "../../../includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbSupersededEntriesCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSupersededEntriesCard::PdbSupersededEntriesCard() {}
PdbSupersededEntriesCard::PdbSupersededEntriesCard(std::string &line)
{
    record_name_ = line.substr(0, 6);
    gmml::Trim(record_name_);

    superseded_date_ = line.substr(11,9);
    gmml::Trim(superseded_date_);

    superseded_id_  = line.substr(31, 48);
    gmml::Trim(superseded_id_);

}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string PdbSupersededEntriesCard::GetRecordName(){
    return record_name_;
}

std::string PdbSupersededEntriesCard::GetSupersededDate(){
    return superseded_date_;
}

std::string PdbSupersededEntriesCard::GetSupersededID(){
    return superseded_id_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSupersededEntriesCard::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

void PdbSupersededEntriesCard::SetSupersededDate(const std::string superseded_date){
    superseded_date_ = superseded_date;
}

void PdbSupersededEntriesCard::SetSupersededID(const std::string superseded_id){
    superseded_id_ = superseded_id;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbSupersededEntriesCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "Superseded Date:" << superseded_date_;
    out << "Superseded ID:" << superseded_id_;
    out << std::endl;
}
