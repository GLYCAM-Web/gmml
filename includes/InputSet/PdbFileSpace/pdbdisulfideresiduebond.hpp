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
            /*! \typedef
              * List of disulfide residues
              */
            typedef std::vector<PdbDisulfideResidue*> DisulfideResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbDisulfideResidueBond();
            /*! \fn
              * Constructor with required parameters
              * @param serial_number Serial number of the disulfide bond
              * @param residues Residues involving in the bond
              * @param bond_length Length of the bond
              */
            PdbDisulfideResidueBond(int serial_number, const DisulfideResidueVector residues, double bond_length);
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbDisulfideResidueBond(std::string& line);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the serial number in a disulfide residue bond
              * @return serial_number_ attribute of the current object of this class
              */
            int GetSerialNumber();
            /*! \fn
              * An accessor function in order to access to the residues in a disulfide residue bond
              * @return residues_ attribute of the current object of this class
              */
            DisulfideResidueVector GetResidues();
            /*! \fn
              * An accessor function in order to access to the bond length in a disulfide residue bond
              * @return bond_length_ attribute of the current object of this class
              */
            double GetBondLength();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the serial number of the current object
              * Set the serial_number_ attribute of the current disulfide residue bond
              * @param serial_number The serial number of the current object
              */
            void SetSerialNumber(int serial_number);
            /*! \fn
              * A mutator function in order to set the residues of the current object
              * Set the residues_ attribute of the current disulfide residue bond
              * @param residues The residues of the current object
              */
            void SetResidues(const DisulfideResidueVector residues);
            /*! \fn
              * A function in order to add the residue to the current object
              * Set the residue_ attribute of the current disulfide residue bond
              * @param residue The residue of the current object
              */
            void AddResidue(PdbDisulfideResidue* residue);
            /*! \fn
              * A mutator function in order to set the bond length of the current object
              * Set the bond_length_ attribute of the current disulfide residue bond
              * @param bond_length The bond length of the current object
              */
            void SetBondLength(double bond_length);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the disulfide residue bond contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            int serial_number_;                 /*!< Serial number of the disulfide bond >*/
            DisulfideResidueVector residues_;   /*!< Residues involving in the disulfide bond >*/
            double bond_length_;                /*!< Length of the bond >*/
    };
}

#endif // PDBDISULFIDERESIDUEBOND_HPP
