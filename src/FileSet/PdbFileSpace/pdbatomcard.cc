
#include "../../../includes/FileSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAtomCard::PdbAtomCard() {}

PdbAtomCard::PdbAtomCard(stringstream &stream_block, string index)
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
        atom->SetAtomCardIndexInResidueSet(index);
        atoms_[atom->GetAtomSerialNumber()] = atom;

        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbAtomCard::GetRecordName()
{
    return record_name_;
}

PdbAtomCard::PdbAtomMap PdbAtomCard::GetAtoms()
{
    return atoms_;
}



//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbAtomCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbAtomCard::SetAtoms(PdbAtomMap atoms)
{
    atoms_.clear();
    for(PdbAtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PdbAtom* atom = (*it).second;
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
void PdbAtomCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "_________________ Atoms _______________" << endl;
    for(PdbAtomCard::PdbAtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
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
