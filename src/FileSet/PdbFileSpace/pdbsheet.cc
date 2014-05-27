// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbsheet.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsheetstrand.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

#include <iostream>

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSheet::PdbSheet() :  sheet_id_(""), number_of_strands_(dNotSet) {}

PdbSheet::PdbSheet(const string &sheet_id, int number_of_strands, const SheetStrandVector strands)
    : sheet_id_(sheet_id), number_of_strands_(number_of_strands)
{
    strands_.clear();
    for(SheetStrandVector::const_iterator it = strands.begin(); it != strands.end(); it++)
    {
        strands_.push_back(*it);
    }
}

PdbSheet::PdbSheet(stringstream& stream_block)
{
    string line;
    bool is_sheet_id_set = false, is_number_of_strands_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_sheet_id_set)
        {
            sheet_id_ = line.substr(11, 3);
            Trim(sheet_id_);
            is_sheet_id_set = true;
        }
        if(!is_number_of_strands_set)
        {
            if(line.substr(14, 2).compare("  ") == 0)
                number_of_strands_ = iNotSet;
            else
                number_of_strands_ = ConvertString<int>(line.substr(14,2));
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
string PdbSheet::GetSheetId()
{
    return sheet_id_;
}

int PdbSheet::GetNumberOfStrands()
{
    return number_of_strands_;
}

PdbSheet::SheetStrandVector PdbSheet::GetStrands()
{
    return strands_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbSheet::SetSheetId(const string sheet_id)
{
    sheet_id_ = sheet_id;
}

void PdbSheet::SetNumberOfStrands(int number_of_strands)
{
    number_of_strands_ = number_of_strands;
}

void PdbSheet::SetStrands(const SheetStrandVector strands)
{
    strands_.clear();
    for(SheetStrandVector::const_iterator it = strands.begin(); it != strands.end(); it++)
    {
        strands_.push_back(*it);
    }
}

void PdbSheet::AddStrand(PdbSheetStrand *strand)
{
    strands_.push_back(strand);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSheet::Print(ostream &out)
{
    out << "Sheet ID: " << sheet_id_
        << ", Number of Strands: ";
    if(number_of_strands_ != iNotSet)
        out << number_of_strands_;
    else
        out << " ";
    out << endl << "------------------- Strands ------------------" << endl;
    for(PdbSheet::SheetStrandVector::iterator it = strands_.begin(); it != strands_.end(); it++)
    {
        (*it)->Print(out);
        out << endl;
    }
    out << endl;
}
