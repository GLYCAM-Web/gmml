// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbformulacard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbformula.hpp"
#include "../../../includes/utils.hpp"


using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFormulaCard::PdbFormulaCard() : record_name_("FORMUL") {}

PdbFormulaCard::PdbFormulaCard(const string &record_name) : record_name_(record_name) {}

PdbFormulaCard::PdbFormulaCard(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }
        stringstream formula_block;
        formula_block << line << endl;
        string heterogen_identifier = line.substr(13,3);

        getline(stream_block, line);

        while (!Trim(line).empty() && line.substr(13,3) == heterogen_identifier){
            formula_block << line << endl;
            getline(stream_block, line);
        }
        PdbFormula* formula = new PdbFormula(formula_block);
        formulas_[heterogen_identifier] = formula;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbFormulaCard::GetRecordName()
{
    return record_name_;
}

PdbFormulaCard::FormulaMap PdbFormulaCard::GetFormulas()
{
    return formulas_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbFormulaCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

