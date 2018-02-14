#include "../../../includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSupersededEntriesCard::PdbSupersededEntriesCard() {}
PdbSupersededEntriesCard::PdbSupersededEntriesCard(string &line)
{
    record_name_ = line.substr(0, 6);
    Trim(record_name_);

    superseded_date_ = line.substr(11,9);
    Trim(superseded_date_);

    superseded_id_  = line.substr(31, 48);
    Trim(superseded_id_);

}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbSupersededEntriesCard::GetRecordName(){
    return record_name_;
}

string PdbSupersededEntriesCard::GetSupersededDate(){
    return superseded_date_;
}

string PdbSupersededEntriesCard::GetSupersededID(){
    return superseded_id_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSupersededEntriesCard::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbSupersededEntriesCard::SetSupersededDate(const string superseded_date){
    superseded_date_ = superseded_date;
}

void PdbSupersededEntriesCard::SetSupersededID(const string superseded_id){
    superseded_id_ = superseded_id;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbSupersededEntriesCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "Superseded Date:" << superseded_date_;
    out << "Superseded ID:" << superseded_id_;
    out << endl;
}
