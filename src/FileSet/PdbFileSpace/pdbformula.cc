// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbformula.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFormula::PdbFormula() : heterogen_identifier_(""), component_number_(kNotSet), chemical_formula_("") {}
PdbFormula::PdbFormula(const string &heterogen_identifier, int component_number, const string &chemical_formula)
    : heterogen_identifier_(heterogen_identifier), component_number_(component_number), chemical_formula_(chemical_formula) {}

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

