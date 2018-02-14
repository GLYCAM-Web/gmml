#include "../../../includes/InputSet/PdbFileSpace/pdbrevisiondatacard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbRevisionDataCard::PdbRevisionDataCard() {}
PdbRevisionDataCard::PdbRevisionDataCard(string &line)
{
    record_name_ = line.substr(0, 6);
    Trim(record_name_);

    string num = line.substr(7,3);
    Trim(num);
    mod_num_ = ConvertString<int>(num);

    mod_date_ = line.substr(13,9);
    Trim(mod_date_);

    mod_id_  = line.substr(23, 4);
    Trim(mod_id_);

    string type = line.substr(32,1);
    mod_type_ = ConvertString<int>(type);

    mod_detail_ = line.substr(39, 40);
    Trim(mod_detail_);
    
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbRevisionDataCard::GetRecordName(){
    return record_name_;
}

int PdbRevisionDataCard::GetModificationNumber(){
    return mod_num_;
}

string PdbRevisionDataCard::GetModificationDate(){
    return mod_date_;
}

string PdbRevisionDataCard::GetModificationID(){
    return mod_id_;
}

int PdbRevisionDataCard::GetModificationType(){
    return mod_type_;
}

string PdbRevisionDataCard::GetModificationDetails(){
    return mod_detail_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbRevisionDataCard::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbRevisionDataCard::SetModificationNumber(const int mod_num){
    mod_num_ = mod_num;
}

void PdbRevisionDataCard::SetModificationDate(const string mod_date){
    mod_date_ = mod_date;
}

void PdbRevisionDataCard::SetModificationID(const string mod_id){
    mod_id_ = mod_id;
}

void PdbRevisionDataCard::SetModificationType(const int mod_type){
    mod_type_ = mod_type;
}

void PdbRevisionDataCard::SetModificationDetails(const string mod_detail){
    mod_detail_ = mod_detail;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbRevisionDataCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "Modification Number:" << mod_num_;
    out << "Modification Date:" << mod_date_;
    out << "Modification ID:" << mod_id_;
    out << "Modification Type:" << mod_type_;
    out << "Modification Details:" << mod_detail_;
    out << endl;
}
