
#include "../../../includes/InputSet/PdbFileSpace/pdbformulacard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbformulasection.hpp"
#include "../../../includes/utils.hpp"


using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFormulaSection::PdbFormulaSection() : record_name_("FORMUL") {}

PdbFormulaSection::PdbFormulaSection(const string &record_name) : record_name_(record_name) {}

PdbFormulaSection::PdbFormulaSection(stringstream &stream_block)
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
        PdbFormulaCard* formula_card = new PdbFormulaCard(formula_block);
        heterogen_identifier = Trim(heterogen_identifier);
        formula_cards_[heterogen_identifier] = formula_card;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbFormulaSection::GetRecordName()
{
    return record_name_;
}

PdbFormulaSection::FormulaCardMap PdbFormulaSection::GetFormulaCards()
{
    return formula_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbFormulaSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbFormulaSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "=========== Formulas =============" << endl;
    for(PdbFormulaSection::FormulaCardMap::iterator it = formula_cards_.begin(); it != formula_cards_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
