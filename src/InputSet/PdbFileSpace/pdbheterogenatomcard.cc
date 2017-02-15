#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenAtomCard::PdbHeterogenAtomCard() : record_name_("HETATM") {}

PdbHeterogenAtomCard::PdbHeterogenAtomCard(stringstream &stream_block, string index)
{
    heterogen_atoms_ = PdbHeterogenAtomMap();
    ordered_heterogen_atoms_ = PdbHeterogenAtomOrderVector();
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
//            record_name_ = line.substr(0,6);
            record_name_ = "HETATM";
            Trim(record_name_);
            is_record_name_set=true;
        }

        PdbAtom* atom = new PdbAtom(line);
        int ch = 97 + ConvertString<int>(Split(index,"_")[1]);
        atom->SetAtomCardIndexInResidueSet(index);
        atom->SetAtomChainId((char)ch);
        heterogen_atoms_[atom->GetAtomSerialNumber()] = atom;
        ordered_heterogen_atoms_.push_back(atom);

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

PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector PdbHeterogenAtomCard::GetOrderedHeterogenAtoms()
{
    return ordered_heterogen_atoms_;
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

void PdbHeterogenAtomCard::SetOrderedHeterogenAtoms(PdbHeterogenAtomOrderVector ordered_heterogen_atoms)
{
    ordered_heterogen_atoms_.clear();
    for(PdbHeterogenAtomOrderVector::iterator it = ordered_heterogen_atoms.begin(); it != ordered_heterogen_atoms.end(); it++)
    {
        PdbAtom* atom = (*it);
        ordered_heterogen_atoms_.push_back(atom);
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
