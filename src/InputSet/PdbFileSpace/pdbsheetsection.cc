// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbsheetcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetsection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSheetSection::PdbSheetSection() : record_name_("SHEET") {}

PdbSheetSection::PdbSheetSection(const string &record_name) : record_name_(record_name) {}

PdbSheetSection::PdbSheetSection(stringstream &stream_block)
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
        stringstream sheet_block;
        sheet_block << line << endl;
        string sheet_id = line.substr(11,3);

        getline(stream_block, line);
        temp = line;

        while (!Trim(temp).empty() && line.substr(11,3) == sheet_id){
            sheet_block << line << endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbSheetCard* sheet = new PdbSheetCard(sheet_block);
        sheet_id = Trim(sheet_id);
        sheet_cards_[sheet_id] = sheet;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbSheetSection::GetRecordName()
{
    return record_name_;
}

PdbSheetSection::SheetCardMap PdbSheetSection::GetSheets()
{
    return sheet_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbSheetSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSheetSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "============== Sheets ==============" << endl;
    for(PdbSheetSection::SheetCardMap::iterator it = sheet_cards_.begin(); it != sheet_cards_.end(); it++)
    {
        out << "Sheet ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
