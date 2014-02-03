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
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }
        stringstream formula_block;
        formula_block << line << endl;
        string heterogen_identifier = line.substr(12,3);

        getline(stream_block, line);
        temp = line;

        while (!Trim(temp).empty() && line.substr(12,3) == heterogen_identifier){
            formula_block << line << endl;
            getline(stream_block, line);
            temp = line;
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
void PdbFormulaCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "=========== Formulas =============" << endl;
    for(PdbFormulaCard::FormulaMap::iterator it = formulas_.begin(); it != formulas_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
