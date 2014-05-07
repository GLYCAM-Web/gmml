#include "../../../includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenAtomCard::PdbHeterogenAtomCard() {}

PdbHeterogenAtomCard::PdbHeterogenAtomCard(stringstream &stream_block)
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

        PdbAtom* atom = new PdbAtom(line);
        heterogen_atoms_[atom->GetAtomSerialNumber()] = atom;

        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHeterogenAtomCard::GetRecordName()
{
    return record_name_;
}

PdbHeterogenAtomCard::PdbHeterogenAtomMap PdbHeterogenAtomCard::GetHeterogenAtoms()
{
    return heterogen_atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenAtomCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbHeterogenAtomCard::SetHeterogenAtoms(PdbHeterogenAtomMap heterogen_atoms)
{
    heterogen_atoms_.clear();
    for(PdbHeterogenAtomMap::iterator it = heterogen_atoms.begin(); it != heterogen_atoms.end(); it++)
    {
        PdbHeterogenAtom* heterogen_atom = (*it).second;
        int serial_number = (*it).first;
        heterogen_atoms_[serial_number] = heterogen_atom;
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenAtomCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "________________ Heterogen Atoms ___________________" << endl;
    for(PdbHeterogenAtomCard::PdbHeterogenAtomMap::iterator it = heterogen_atoms_.begin(); it != heterogen_atoms_.end(); it++)
    {
        out << "Atom Serial Number: ";
        if((it)->first != iNotSet)
            out << (it)->first << endl;
        else
            out << " " << endl;
        (it)->second->Print();
        out << endl;
    }
}
