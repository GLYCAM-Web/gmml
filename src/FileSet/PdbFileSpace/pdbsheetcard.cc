// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbsheetcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsheet.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSheetCard::PdbSheetCard() : record_name_("SHEET") {}

PdbSheetCard::PdbSheetCard(const string &record_name) : record_name_(record_name) {}

PdbSheetCard::PdbSheetCard(stringstream &stream_block)
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
        PdbSheet* sheet = new PdbSheet(sheet_block);
        sheet_id = Trim(sheet_id);
        sheets_[sheet_id] = sheet;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbSheetCard::GetRecordName()
{
    return record_name_;
}

PdbSheetCard::SheetMap PdbSheetCard::GetSheets()
{
    return sheets_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbSheetCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSheetCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "============== Sheets ==============" << endl;
    for(PdbSheetCard::SheetMap::iterator it = sheets_.begin(); it != sheets_.end(); it++)
    {
        out << "Sheet ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
