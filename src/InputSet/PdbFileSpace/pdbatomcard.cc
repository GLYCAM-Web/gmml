
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAtomCard::PdbAtomCard() : record_name_("ATOM") {}

PdbAtomCard::PdbAtomCard(stringstream &stream_block, string index)
{
    atoms_ = PdbAtomMap();
    ordered_atoms_ = PdbAtomOrderVector();
    string line;
    bool is_record_name_set = false;
//    cout << stream_block.str() << endl;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
//            record_name_ = line.substr(0,6);
            record_name_ = "ATOM";
            Trim(record_name_);
            is_record_name_set=true;
        }

        PdbAtom* atom = new PdbAtom(line);
        int ch = 65 + ConvertString<int>(Split(index, "_")[1]);
        atom->SetAtomCardIndexInResidueSet(index);
        atom->SetAtomChainId((char)ch);
        atoms_[atom->GetAtomSerialNumber()] = atom;
        ordered_atoms_.push_back(atom);

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
PdbAtomCard::PdbAtomOrderVector PdbAtomCard::GetOrderedAtoms()
{
    return ordered_atoms_;
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
void PdbAtomCard::SetOrderedAtoms(PdbAtomOrderVector ordered_atoms)
{
    ordered_atoms_.clear();
    for(PdbAtomOrderVector::iterator it = ordered_atoms.begin(); it != ordered_atoms.end(); it++)
    {
        PdbAtom* atom = (*it);
        ordered_atoms_.push_back(atom);
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
