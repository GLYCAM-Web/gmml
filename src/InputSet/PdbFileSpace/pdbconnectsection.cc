#include "../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbConnectSection::PdbConnectSection() : record_name_("CONECT")
{
    bonded_atom_serial_numbers_ = BondedAtomsSerialNumbersMap();
}

PdbConnectSection::PdbConnectSection(stringstream &stream_block)
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

        int atom_serial_number;
        if(line.substr(6, 5) != "     ")
            atom_serial_number = ConvertString<int>(line.substr(6,5));
        else
            atom_serial_number = iNotSet;
        vector<int> bonded_atom_serial_numbers;
        if(line.substr(11,5) != "     ")
            bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(11,5)));
        if(line.substr(16,5) != "     ")
            bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(16,5)));
        if(line.substr(21,5) != "     ")
            bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(21,5)));
        if(line.substr(26,5) != "     ")
            bonded_atom_serial_numbers.push_back(ConvertString<int>(line.substr(26,5)));
        BondedAtomsSerialNumbersMap::iterator it = bonded_atom_serial_numbers_.find(atom_serial_number);
        if(it->first == atom_serial_number)
        {
            for(vector<int>::iterator it = bonded_atom_serial_numbers.begin(); it != bonded_atom_serial_numbers.end(); it++)
                bonded_atom_serial_numbers_.find(atom_serial_number)->second.push_back((*it));
        }
        else
        {
            bonded_atom_serial_numbers_[atom_serial_number] = bonded_atom_serial_numbers;
        }
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
string PdbConnectSection::GetRecordName()
{
    return record_name_;
}

PdbConnectSection::BondedAtomsSerialNumbersMap PdbConnectSection::GetBondedAtomsSerialNumbers()
{
    return bonded_atom_serial_numbers_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbConnectSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbConnectSection::SetBondedAtomsSerialNumbers(BondedAtomsSerialNumbersMap bonded_atom_serial_numbers)
{
    bonded_atom_serial_numbers_.clear();
    for(BondedAtomsSerialNumbersMap::iterator it = bonded_atom_serial_numbers.begin(); it != bonded_atom_serial_numbers.end(); it++)
    {
        int source_serial_number = (*it).first;
        vector<int> bonded_serial_numbers = (*it).second;
        bonded_atom_serial_numbers_[source_serial_number] = bonded_serial_numbers;
    }
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbConnectSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "================= Bonded Atoms Serial Numbers ================" << endl;
    for(PdbConnectSection::BondedAtomsSerialNumbersMap::iterator it = bonded_atom_serial_numbers_.begin(); it != bonded_atom_serial_numbers_.end(); it++)
    {
        out << "Atom Serial Number: ";
        if((it)->first != iNotSet)
            out << (it)->first << endl;
        else
            out << " " << endl;
        for(unsigned int i = 0; i < (it)->second.size(); i++)
        {
            out << (it)->second.at(i) << ", ";
        }
        out << endl;
    }
    out << endl;
}
