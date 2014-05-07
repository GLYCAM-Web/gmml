#ifndef PDBPREPROCESSORREPLACEDHYDROGEN_HPP
#define PDBPREPROCESSORREPLACEDHYDROGEN_HPP

#include <string>
#include <iostream>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorReplacedHydrogen
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorReplacedHydrogen();
            PdbPreprocessorReplacedHydrogen(char residue_chain_id, int atom_serial_number, std::string atom_name, std::string residue_name, int residue_sequence_number);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the atom serial number
              * @return atom_serial_number_ attribute of the current object of this class
              */
            int GetAtomSerialNumber();
            /*! \fn
              * An accessor function in order to access to the atom name
              * @return atom_name_ attribute of the current object of this class
              */
            std::string GetAtomName();
            /*! \fn
              * An accessor function in order to access to the residue name
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the residue sequence number
              * @return residue_sequence_number_ attribute of the current object of this class
              */
            int GetResidueSequenceNumber();
            char GetResidueChainId();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom serial number of the current object
              * Set the atom_serial_number_ attribute of the current pdb preprocessor replaced hydrogen
              * @param atom_seial_number The atom serial number attribute of the current object
              */
            void SetAtomSerialNumber(int atom_serial_number);
            /*! \fn
              * A mutator function in order to set the atom name of the current object
              * Set the atom_name_ attribute of the current pdb preprocessor replaced hydrogen
              * @param atom_name The atom name attribute of the current object
              */
            void SetAtomName(std::string atom_name);
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current pdb preprocessor replaced hydrogen
              * @param residue_name The residue name attribute of the current object
              */
            void SetResidueName(std::string residue_name);
            /*! \fn
              * A mutator function in order to set the residue sequence number of the current object
              * Set the residue_sequence_number_ attribute of the current pdb preprocessor replaced hydrogen
              * @param residue_sequence_number The residue sequence number attribute of the current object
              */
            void SetResidueSequenceNumber(int residue_sequence_number);
            void SetResidueChainId(char residue_chain_id);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor replaced hydrogen contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int atom_serial_number_;
            std::string atom_name_;
            std::string residue_name_;
            int residue_sequence_number_;
            char residue_chain_id_;

    };
}


#endif // PDBPREPROCESSORREPLACEDHYDROGEN_HPP
