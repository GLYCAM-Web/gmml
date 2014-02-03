// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBDISULFIDERESIDUEBOND_HPP
#define PDBDISULFIDERESIDUEBOND_HPP

#include <vector>
#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbDisulfideResidue;
    class PdbDisulfideResidueBond
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<PdbDisulfideResidue*> DisulfideResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbDisulfideResidueBond();
            PdbDisulfideResidueBond(int serial_number, const DisulfideResidueVector residues, double bond_length);
            PdbDisulfideResidueBond(std::string& line);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            int GetSerialNumber();
            DisulfideResidueVector GetResidues();
            double GetBondLength();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetSerialNumber(int serial_number);
            void SetResidues(const DisulfideResidueVector residues);
            void AddResidue(PdbDisulfideResidue* residue);
            void SetBondLength(double bond_length);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            int serial_number_;
            DisulfideResidueVector residues_;
            double bond_length_;
    };
}

#endif // PDBDISULFIDERESIDUEBOND_HPP
