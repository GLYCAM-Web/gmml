// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDisulfideResidueBond::PdbDisulfideResidueBond() : serial_number_(dNotSet), bond_length_(dNotSet) {}

PdbDisulfideResidueBond::PdbDisulfideResidueBond(int serial_number, const DisulfideResidueVector residues, double bond_length)
    : serial_number_(serial_number), bond_length_(bond_length)
{
    residues_.clear();
    for(DisulfideResidueVector::const_iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}

PdbDisulfideResidueBond::PdbDisulfideResidueBond(string &line)
{
    string temp0;
    char temp1, temp2;
    int temp3, temp4;

    temp0 = line.substr(11,3);
    temp0 = Trim(temp0);
    if(line.substr(15,1) == " ")
        temp1 = ' ';
    else
        temp1 = ConvertString<char>(line.substr(15, 1));
    if(line.substr(21,1) == " ")
        temp2 = ' ';
    else
        temp2 = ConvertString<char>(line.substr(21, 1));
    if(line.substr(17, 4) == "    ")
        temp3 = iNotSet;
    else
        temp3 = ConvertString<int>(line.substr(17, 4));
    if(line.substr(59,6) == "      ")
        temp4 = iNotSet;
    else
        temp4 = ConvertString<int>(line.substr(59,6));
    PdbDisulfideResidue* residue_1 = new PdbDisulfideResidue(temp0, temp1, temp3, temp2, temp4);

    temp0 = line.substr(25, 3);
    temp0 = Trim(temp0);
    if(line.substr(29,1) == " ")
        temp1 = ' ';
    else
        temp1 = ConvertString<char>(line.substr(29, 1));
    if(line.substr(35,1) == " ")
        temp2 = ' ';
    else
        temp2 = ConvertString<char>(line.substr(35, 1));
    if(line.substr(31, 4) == "    ")
        temp3 = iNotSet;
    else
        temp3 = ConvertString<int>(line.substr(31, 4));
    if(line.substr(66, 6) == "      ")
        temp4 = iNotSet;
    else
        temp4 = ConvertString<int>(line.substr(66,6));

    PdbDisulfideResidue* residue_2 = new PdbDisulfideResidue(temp0, temp1, temp3, temp2, temp4);

    residues_.push_back(residue_1);
    residues_.push_back(residue_2);

    if(line.substr(7, 3) == "   ")
        serial_number_ = iNotSet;
    else
        serial_number_ = ConvertString<int>(line.substr(7,3));
    if(line.substr(73, 5) == "     ")
        bond_length_ = dNotSet;
    else
        bond_length_ = ConvertString<double>(line.substr(73,5));
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int PdbDisulfideResidueBond::GetSerialNumber()
{
    return serial_number_;
}

PdbDisulfideResidueBond::DisulfideResidueVector PdbDisulfideResidueBond::GetResidues()
{
    return residues_;
}

double PdbDisulfideResidueBond::GetBondLength()
{
    return bond_length_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbDisulfideResidueBond::SetSerialNumber(int serial_number)
{
    serial_number_ = serial_number;
}

void PdbDisulfideResidueBond::SetResidues(const DisulfideResidueVector residues)
{
    residues_.clear();
    for(DisulfideResidueVector::const_iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}

void PdbDisulfideResidueBond::AddResidue(PdbDisulfideResidue* residue)
{
    residues_.push_back(residue);
}

void PdbDisulfideResidueBond::SetBondLength(double bond_length)
{
    bond_length_ = bond_length;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbDisulfideResidueBond::Print(ostream &out)
{
    out << "Serial Number: ";
    if(serial_number_ != iNotSet)
        out << serial_number_ << endl;
    else
        out << " " << endl;
    out << "----------------- Residues -----------------" << endl;
    for(PdbDisulfideResidueBond::DisulfideResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        (*it)->Print(out);
        out << endl;
    }
    out << "Bond Length: ";
    if(bond_length_ != dNotSet)
        out << bond_length_ << endl << endl;
    else
        out << " " << endl << endl;
}
