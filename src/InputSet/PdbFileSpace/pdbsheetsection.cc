// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbsheetcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbSheetSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSheetSection::PdbSheetSection() : record_name_("SHEET") {}

PdbSheetSection::PdbSheetSection(const std::string &record_name) : record_name_(record_name) {}

PdbSheetSection::PdbSheetSection(std::stringstream &stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
        std::stringstream sheet_block;
        sheet_block << line << std::endl;
        std::string sheet_id = line.substr(11,3);

        getline(stream_block, line);
        temp = line;

        while (!gmml::Trim(temp).empty() && line.substr(11,3) == sheet_id){
            sheet_block << line << std::endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbSheetCard* sheet = new PdbSheetCard(sheet_block);
        sheet_id = gmml::Trim(sheet_id);
        sheet_cards_[sheet_id] = sheet;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbSheetSection::GetRecordName()
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
void PdbSheetSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSheetSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "============== Sheets ==============" << std::endl;
    for(PdbSheetSection::SheetCardMap::iterator it = sheet_cards_.begin(); it != sheet_cards_.end(); it++)
    {
        out << "Sheet ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
