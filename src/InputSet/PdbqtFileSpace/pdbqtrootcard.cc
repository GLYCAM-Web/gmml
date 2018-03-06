#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtRootCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtRootCard::PdbqtRootCard() : record_name_("ROOT")
{
    root_atoms_ = NULL;
}

PdbqtRootCard::PdbqtRootCard(std::stringstream &root_block)
{
    root_atoms_ = NULL;
    std::string line;
    bool is_record_name_set = false;
    getline(root_block, line);
    std::string temp = line;
    if(line.find("ROOT") != std::string::npos)
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
    }

    getline(root_block,line);
    std::stringstream stream_block;
    while(line.find("ATOM") != std::string::npos || line.find("HETATOM") != std::string::npos)
    {
        stream_block << line << std::endl;
        getline(root_block,line);
        temp = line;
    }
    if(line.find("ENDROOT") != std::string::npos)
    {
        root_atoms_ = new PdbqtFileSpace::PdbqtAtomCard(stream_block);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbqtRootCard::GetRecordName()
{
    return record_name_;
}
PdbqtFileSpace::PdbqtAtomCard* PdbqtRootCard::GetRootAtoms()
{
    return root_atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtRootCard::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbqtRootCard::SetRootAtoms(PdbqtFileSpace::PdbqtAtomCard* root_atoms)
{
    root_atoms_ = new PdbqtFileSpace::PdbqtAtomCard();
    root_atoms_->SetRecordName(root_atoms->GetRecordName());
    root_atoms_->SetAtoms(root_atoms->GetAtoms());
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtRootCard::Print(std::ostream &out)
{
    if(root_atoms_ != NULL)
        root_atoms_->Print(out);
}
