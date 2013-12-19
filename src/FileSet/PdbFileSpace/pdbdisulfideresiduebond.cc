// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/common.hpp"

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

