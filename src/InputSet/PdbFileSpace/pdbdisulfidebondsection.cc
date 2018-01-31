
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfidebondsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDisulfideBondSection::PdbDisulfideBondSection() : record_name_("SSBOND") {}

PdbDisulfideBondSection::PdbDisulfideBondSection(const string &record_name) : record_name_(record_name) {}

PdbDisulfideBondSection::PdbDisulfideBondSection(stringstream &stream_block)
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
string PdbDisulfideBondSection::GetRecordName()
{
    return record_name_;
}

PdbDisulfideBondSection::DisulfideResidueBondMap PdbDisulfideBondSection::GetDisulfideResidueBonds()
{
    return disulfide_residue_bonds_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbDisulfideBondSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbDisulfideBondSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "================= Disulfide Bond =================" << endl;
    for(PdbDisulfideBondSection::DisulfideResidueBondMap::iterator it = disulfide_residue_bonds_.begin(); it != disulfide_residue_bonds_.end(); it++)
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
