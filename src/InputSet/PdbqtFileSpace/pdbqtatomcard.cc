#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtAtomCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtAtomCard::PdbqtAtomCard() : record_name_("ATOM"){}

PdbqtAtomCard::PdbqtAtomCard(std::stringstream &stream_block)
{
    atoms_ = PdbqtAtomMap();
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = "ATOM";
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }

        PdbqtAtom* atom = new PdbqtAtom(line);
        atoms_[atom->GetAtomSerialNumber()] = atom;

        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbqtAtomCard::GetRecordName()
{
    return record_name_;
}

PdbqtAtomCard::PdbqtAtomMap PdbqtAtomCard::GetAtoms()
{
    return atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtAtomCard::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbqtAtomCard::SetAtoms(PdbqtAtomMap atoms)
{
    atoms_.clear();
    for(PdbqtAtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PdbqtAtom* atom = (*it).second;
        int serial_number = (*it).first;
        atoms_[serial_number] = atom;
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtAtomCard::Print(std::ostream &out)
{
//    out << "Record Name: " << record_name_ << std::endl <<
    out << "_________________ Atoms _______________" << std::endl;
    for(PdbqtAtomCard::PdbqtAtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        out << "Atom Serial Number: ";
        if((it)->first != gmml::iNotSet)
            out << (it)->first << std::endl;
        else
            out << " " << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
    out << std::endl;
}
