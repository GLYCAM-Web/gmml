// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbsheetcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

#include <iostream>

using PdbFileSpace::PdbSheetCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSheetCard::PdbSheetCard() :  sheet_id_(""), number_of_strands_(gmml::dNotSet) {}

PdbSheetCard::PdbSheetCard(const std::string &sheet_id, int number_of_strands, const SheetStrandVector strands)
    : sheet_id_(sheet_id), number_of_strands_(number_of_strands)
{
    strands_.clear();
    for(SheetStrandVector::const_iterator it = strands.begin(); it != strands.end(); it++)
    {
        strands_.push_back(*it);
    }
}

PdbSheetCard::PdbSheetCard(std::stringstream& stream_block)
{
    std::string line;
    bool is_sheet_id_set = false, is_number_of_strands_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_sheet_id_set)
        {
            sheet_id_ = line.substr(11, 3);
            gmml::Trim(sheet_id_);
            is_sheet_id_set = true;
        }
        if(!is_number_of_strands_set)
        {
            if(line.substr(14, 2).compare("  ") == 0)
                number_of_strands_ = gmml::iNotSet;
            else
                number_of_strands_ = gmml::ConvertString<int>(line.substr(14,2));
            is_number_of_strands_set = true;
        }

        PdbSheetStrand* sheet_strand = new PdbSheetStrand(line);

        strands_.push_back(sheet_strand);

        getline(stream_block, line);
        temp = line;
    }
}


//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbSheetCard::GetSheetId()
{
    return sheet_id_;
}

int PdbSheetCard::GetNumberOfStrands()
{
    return number_of_strands_;
}

PdbSheetCard::SheetStrandVector PdbSheetCard::GetStrands()
{
    return strands_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbSheetCard::SetSheetId(const std::string sheet_id)
{
    sheet_id_ = sheet_id;
}

void PdbSheetCard::SetNumberOfStrands(int number_of_strands)
{
    number_of_strands_ = number_of_strands;
}

void PdbSheetCard::SetStrands(const SheetStrandVector strands)
{
    strands_.clear();
    for(SheetStrandVector::const_iterator it = strands.begin(); it != strands.end(); it++)
    {
        strands_.push_back(*it);
    }
}

void PdbSheetCard::AddStrand(PdbSheetStrand *strand)
{
    strands_.push_back(strand);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSheetCard::Print(std::ostream &out)
{
    out << "Sheet ID: " << sheet_id_
        << ", Number of Strands: ";
    if(number_of_strands_ != gmml::iNotSet)
        out << number_of_strands_;
    else
        out << " ";
    out << std::endl << "------------------- Strands ------------------" << std::endl;
    for(PdbSheetCard::SheetStrandVector::iterator it = strands_.begin(); it != strands_.end(); it++)
    {
        (*it)->Print(out);
        out << std::endl;
    }
    out << std::endl;
}
