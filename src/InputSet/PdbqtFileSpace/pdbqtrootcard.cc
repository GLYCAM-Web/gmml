#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtRootCard::PdbqtRootCard() : record_name_("ROOT")
{
    root_atoms_ = NULL;
}

PdbqtRootCard::PdbqtRootCard(stringstream &root_block)
{
    root_atoms_ = NULL;
    string line;
    bool is_record_name_set = false;
    getline(root_block, line);
    string temp = line;
    if(line.find("ROOT") != string::npos)
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
    }

    getline(root_block,line);
    stringstream stream_block;
    while(line.find("ATOM") != string::npos || line.find("HETATOM") != string::npos)
    {
        stream_block << line << endl;
        getline(root_block,line);
        temp = line;
    }
    if(line.find("ENDROOT") != string::npos)
    {
        root_atoms_ = new PdbqtAtomCard(stream_block);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbqtRootCard::GetRecordName()
{
    return record_name_;
}
PdbqtAtomCard* PdbqtRootCard::GetRootAtoms()
{
    return root_atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtRootCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbqtRootCard::SetRootAtoms(PdbqtAtomCard* root_atoms)
{    
    root_atoms_ = new PdbqtAtomCard();
    root_atoms_->SetRecordName(root_atoms->GetRecordName());
    root_atoms_->SetAtoms(root_atoms->GetAtoms());
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtRootCard::Print(ostream &out)
{
    if(root_atoms_ != NULL)
        root_atoms_->Print(out);
}


