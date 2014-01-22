// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfidebondcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDisulfideBondCard::PdbDisulfideBondCard() : record_name_("SSBOND") {}

PdbDisulfideBondCard::PdbDisulfideBondCard(const string &record_name) : record_name_(record_name) {}

PdbDisulfideBondCard::PdbDisulfideBondCard(stringstream &stream_block)
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

        PdbDisulfideResidueBond* disulfide_bond = new PdbDisulfideResidueBond(line);
        disulfide_residue_bonds_[ConvertString<int>(line.substr(7,3))] = disulfide_bond;
        getline(stream_block, line);
    }
}

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

