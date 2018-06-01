
#include "../../../includes/InputSet/PdbFileSpace/pdbformulacard.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbFormulaCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFormulaCard::PdbFormulaCard() : heterogen_identifier_(""), component_number_(gmml::dNotSet), chemical_formula_("") {}
PdbFormulaCard::PdbFormulaCard(const std::string &heterogen_identifier, int component_number, const std::string &chemical_formula)
    : heterogen_identifier_(heterogen_identifier), component_number_(component_number), chemical_formula_(chemical_formula) {}

PdbFormulaCard::PdbFormulaCard(std::stringstream& stream_block)
{
    std::string line;
    bool is_heterogen_identifier_set = false, is_component_number_set = false;
    std::stringstream ss;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_heterogen_identifier_set){
            heterogen_identifier_ = line.substr(12,3);
            gmml::Trim(heterogen_identifier_);
            is_heterogen_identifier_set = true;
        }

        if(!is_component_number_set){
            if(line.substr(8, 2) == "  ")
                component_number_ = gmml::iNotSet;
            else
                component_number_ = gmml::ConvertString<int>(line.substr(8,2));
            is_component_number_set = true;
        }

        ss << line.substr(19,51) << " ";

        getline(stream_block, line);
        temp = line;
    }
    chemical_formula_ = ss.str();
    chemical_formula_ = gmml::Trim(chemical_formula_);
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbFormulaCard::GetHeterogenIdentifier()
{
    return heterogen_identifier_;
}

int PdbFormulaCard::GetComponentNumber()
{
    return component_number_;
}

std::string PdbFormulaCard::GetChemicalFormula()
{
    return chemical_formula_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbFormulaCard::SetHeterogenIdentifier(const std::string heterogen_identifier)
{
    heterogen_identifier_ = heterogen_identifier;
}

void PdbFormulaCard::SetComponentNumber(int component_number)
{
    component_number_ = component_number;
}

void PdbFormulaCard::SetChemicalFormula(const std::string chemical_formula)
{
    chemical_formula_ = chemical_formula;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void  PdbFormulaCard::Print(std::ostream &out)
{
    out << "Heterogen ID: " << heterogen_identifier_
        << ", Component Number: ";
    if(component_number_ != gmml::iNotSet)
        out << component_number_;
    else
        out << " ";
    out << "Chemical Formula: " << chemical_formula_ << std::endl;
}
