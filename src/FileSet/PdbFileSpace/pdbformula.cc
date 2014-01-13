// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbformula.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFormula::PdbFormula() : heterogen_identifier_(""), component_number_(kNotSet), chemical_formula_("") {}
PdbFormula::PdbFormula(const string &heterogen_identifier, int component_number, const string &chemical_formula)
    : heterogen_identifier_(heterogen_identifier), component_number_(component_number), chemical_formula_(chemical_formula) {}

PdbFormula::PdbFormula(stringstream& stream_block)
{
    string line;
    bool is_heterogen_identifier_set = false, is_component_number_set = false;
    stringstream ss;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_heterogen_identifier_set){
            heterogen_identifier_ = line.substr(12,3);
            is_heterogen_identifier_set = true;
        }

        if(!is_component_number_set){
            component_number_ = ConvertString<int>(line.substr(8,2));
            is_component_number_set = true;
        }

        ss << line.substr(19,51) << " ";

        getline(stream_block, line);
    }
    chemical_formula_ = ss.str();
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbFormula::GetHeterogenIdentifier()
{
    return heterogen_identifier_;
}

int PdbFormula::GetComponentNumber()
{
    return component_number_;
}

string PdbFormula::GetChemicalFormula()
{
    return chemical_formula_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbFormula::SetHeterogenIdentifier(const string heterogen_identifier)
{
    heterogen_identifier_ = heterogen_identifier;
}

void PdbFormula::SetComponentNumber(int component_number)
{
    component_number_ = component_number;
}

void PdbFormula::SetChemicalFormula(const string chemical_formula)
{
    chemical_formula_ = chemical_formula;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

