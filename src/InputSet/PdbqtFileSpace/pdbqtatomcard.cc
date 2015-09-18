#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtAtomCard::PdbqtAtomCard() : record_name_("ATOM"){}

PdbqtAtomCard::PdbqtAtomCard(stringstream &stream_block)
{
    atoms_ = PdbqtAtomMap();
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = "ATOM";
            Trim(record_name_);
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
string PdbqtAtomCard::GetRecordName()
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
void PdbqtAtomCard::SetRecordName(const string record_name)
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
void PdbqtAtomCard::Print(ostream &out)
{
//    out << "Record Name: " << record_name_ << endl <<
    out << "_________________ Atoms _______________" << endl;
    for(PdbqtAtomCard::PdbqtAtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        out << "Atom Serial Number: ";
        if((it)->first != iNotSet)
            out << (it)->first << endl;
        else
            out << " " << endl;
        (it)->second->Print();
        out << endl;
    }
    out << endl;
}
