// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfidebondcard.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDisulfideBondCard::PdbDisulfideBondCard() : record_name_("SSBOND") {}

PdbDisulfideBondCard::PdbDisulfideBondCard(const string &record_name) : record_name_(record_name) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbDisulfideBondCard::GetRecordName()
{
    return record_name_;
}

PdbDisulfideBondCard::DisulfideResidueBondMap PdbDisulfideBondCard::GetDisulfideResidueBonds()
{
    return disulfide_residue_bonds_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbDisulfideBondCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

