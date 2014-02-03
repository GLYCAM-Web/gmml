#include "../../../includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/utils.hpp"

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
        out << "Atom Serial Number: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
