
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfidebondsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbDisulfideBondSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDisulfideBondSection::PdbDisulfideBondSection() : record_name_("SSBOND") {}

PdbDisulfideBondSection::PdbDisulfideBondSection(const std::string &record_name) : record_name_(record_name) {}

PdbDisulfideBondSection::PdbDisulfideBondSection(std::stringstream &stream_block)
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

        PdbDisulfideResidueBond* disulfide_bond = new PdbDisulfideResidueBond(line);
        disulfide_residue_bonds_[disulfide_bond->GetSerialNumber()] = disulfide_bond;
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbDisulfideBondSection::GetRecordName()
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
void PdbDisulfideBondSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbDisulfideBondSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "================= Disulfide Bond =================" << std::endl;
    for(PdbDisulfideBondSection::DisulfideResidueBondMap::iterator it = disulfide_residue_bonds_.begin(); it != disulfide_residue_bonds_.end(); it++)
    {
        out << "Serial Number: ";
        if((it)->first != gmml::iNotSet)
            out << (it)->first << std::endl;
        else
            out << " " << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
