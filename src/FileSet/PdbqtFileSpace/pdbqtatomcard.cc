#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtAtomCard::PdbqtAtomCard(){}

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
}
