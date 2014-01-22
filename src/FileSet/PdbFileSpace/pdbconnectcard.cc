#include "../../../includes/FileSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbConnectCard::PdbConnectCard() {}

PdbConnectCard::PdbConnectCard(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }

        int atom_serial_number = ConvertString<int>(line.substr(6,5));
        vector<int> bonded_atom_serial_numbers;
        bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(11,5)));
        bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(16,5)));
        bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(21,5)));
        bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(26,5)));
        bonded_atom_serial_numbers_[atom_serial_number] = bonded_atom_serial_numbers;
        getline(stream_block, line);
    }
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbConnectCard::GetRecordName()
{
    return record_name_;
}

PdbConnectCard::BondedAtomsSerialNumbersMap PdbConnectCard::GetBondedAtomsSerialNumbers()
{
    return bonded_atom_serial_numbers_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbConnectCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////

