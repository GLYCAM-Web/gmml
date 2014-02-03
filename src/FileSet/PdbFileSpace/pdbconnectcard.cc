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
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }

        int atom_serial_number = ConvertString<int>(line.substr(6,5));
        vector<int> bonded_atom_serial_numbers;
        if(line.substr(11,5) != "     ")
            bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(11,5)));
        if(line.substr(16,5) != "     ")
            bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(16,5)));
        if(line.substr(21,5) != "     ")
            bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(21,5)));
        if(line.substr(26,5) != "     ")
            bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(26,5)));
        bonded_atom_serial_numbers_[atom_serial_number] = bonded_atom_serial_numbers;
        getline(stream_block, line);
        temp = line;
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
void PdbConnectCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "================= Bonded Atoms Serial Numbers ================" << endl;
    for(PdbConnectCard::BondedAtomsSerialNumbersMap::iterator it = bonded_atom_serial_numbers_.begin(); it != bonded_atom_serial_numbers_.end(); it++)
    {
        out << "Atom Serial Number: " << (it)->first << endl;
        for(unsigned int i = 0; i < (it)->second.size(); i++)
        {
            out << (it)->second.at(i) << ", ";
        }
        out << endl;
    }
    out << endl;
}
