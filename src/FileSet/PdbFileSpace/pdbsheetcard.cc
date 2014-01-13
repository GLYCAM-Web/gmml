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
        PdbSheet* sheet = new PdbSheet(ss);
        sheets_[line.substr(11,3)] = sheet;
        getline(stream_block, line);
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

