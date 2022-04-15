
#include "../../../includes/InputSet/PdbFileSpace/pdbformulacard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbformulasection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbFormulaSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFormulaSection::PdbFormulaSection() : record_name_("FORMUL") {}

PdbFormulaSection::PdbFormulaSection(const std::string &record_name) : record_name_(record_name) {}

PdbFormulaSection::PdbFormulaSection(std::stringstream &stream_block)
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
        std::stringstream formula_block;
        formula_block << line << std::endl;
        std::string heterogen_identifier = line.substr(12,3);

        getline(stream_block, line);
        temp = line;

        while (!gmml::Trim(temp).empty() && line.substr(12,3) == heterogen_identifier){
            formula_block << line << std::endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbFormulaCard* formula_card = new PdbFormulaCard(formula_block);
        heterogen_identifier = gmml::Trim(heterogen_identifier);
        formula_cards_[heterogen_identifier] = formula_card;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbFormulaSection::GetRecordName()
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
void PdbFormulaSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbFormulaSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "=========== Formulas =============" << std::endl;
    for(PdbFormulaSection::FormulaCardMap::iterator it = formula_cards_.begin(); it != formula_cards_.end(); it++)
    {
        out << "Heterogen ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
