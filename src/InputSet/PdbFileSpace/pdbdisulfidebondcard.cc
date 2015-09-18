// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfidebondcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

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
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }

        PdbDisulfideResidueBond* disulfide_bond = new PdbDisulfideResidueBond(line);
        disulfide_residue_bonds_[disulfide_bond->GetSerialNumber()] = disulfide_bond;
        getline(stream_block, line);
        temp = line;
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
void PdbDisulfideBondCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "================= Disulfide Bond =================" << endl;
    for(PdbDisulfideBondCard::DisulfideResidueBondMap::iterator it = disulfide_residue_bonds_.begin(); it != disulfide_residue_bonds_.end(); it++)
    {
        out << "Serial Number: ";
        if((it)->first != iNotSet)
            out << (it)->first << endl;
        else
            out << " " << endl;
        (it)->second->Print();
        out << endl;
    }
}
