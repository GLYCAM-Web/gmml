// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfideresidue.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDisulfideResidueBond::PdbDisulfideResidueBond() : serial_number_(kNotSet), bond_length_(kNotSet) {}

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
    PdbDisulfideResidue* residue_1 = new PdbDisulfideResidue(line.substr(11, 3), ConvertString<char>(line.substr(15, 1)), ConvertString<int>(line.substr(17, 4)), ConvertString<char>(line.substr(21, 1)), ConvertString<int>(line.substr(59,6)) );
    PdbDisulfideResidue* residue_2 = new PdbDisulfideResidue(line.substr(25, 3), ConvertString<char>(line.substr(29, 1)), ConvertString<int>(line.substr(31, 4)), ConvertString<char>(line.substr(35, 1)), ConvertString<int>(line.substr(66,6)));

    residues_.push_back(residue_1);
    residues_.push_back(residue_2);
    serial_number_ = ConvertString<int>(line.substr(7,3));
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

