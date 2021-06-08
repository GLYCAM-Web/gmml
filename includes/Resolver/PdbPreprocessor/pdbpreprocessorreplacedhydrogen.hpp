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
            /*! \fn
              * Constructor in order to initialize the attributes of an entity of removed/replaced hydrogen atom
              * @param residue_chain_id Chain id of the residue that the removed/replaced hydrogen atom belongs to
              * @param atom_serial_number Serial number of the removed/replaced hydrogen atom
              * @param atom_name Name of the removed/replaced hydrogen atom
              * @param residue_name Name of the residue that the removed/replaced hydrogen atom belongs to
              * @param residue_sequence_number Sequence number of the residue that the removed/replaced hydrogen atom belongs to
              * @param residue_insertion_code Insertion code of the residue that the removed/replaced hydrogen atom belongs to
              * @param residue_alternate_location Alternate location of the residue that the removed/replaced hydrogen atom belongs to
              */
            PdbPreprocessorReplacedHydrogen(char residue_chain_id, int atom_serial_number, std::string atom_name, std::string residue_name, int residue_sequence_number,
                                            char residue_insertion_code, char residue_alternate_location);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
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
            /*! \fn
              * An accessor function in order to access to the residue chain id
              * @return residue_chain_id_ attribute of the current object of this class
              */
            char GetResidueChainId();
            /*! \fn
              * An accessor function in order to access to the residue insertion code
              * @return residue_insertion_code_ attribute of the current object of this class
              */
            char GetResidueInsertionCode();
            /*! \fn
              * An accessor function in order to access to the residue alternate location
              * @return residue_alternate_location_ attribute of the current object of this class
              */
            char GetResidueAlternateLocation();
/** @*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
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
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current pdb preprocessor replaced hydrogen
              * @param residue_chain_id The residue chain id attribute of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the residue insertion code of the current object
              * Set the residue_insertion_code_ attribute of the current pdb preprocessor replaced hydrogen
              * @param residue_insertion_code The residue insertion code attribute of the current object
              */
            void SetResidueInsertionCode(char residue_insertion_code);
            /*! \fn
              * A mutator function in order to set the residue alternate location of the current object
              * Set the residue_alternate_location_ attribute of the current pdb preprocessor replaced hydrogen
              * @param residue_alternate_location The residue alternate location attribute of the current object
              */
            void SetResidueAlternateLocation(char residue_alternate_location);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor replaced hydrogen contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int atom_serial_number_;                /*!< Serial number of the removed/replaced hydrogen atom >*/
            std::string atom_name_;                 /*!< Name of the removed/replaced hydrogen atom >*/
            std::string residue_name_;              /*!< Name of the residue of the removed/replaced hydrogen atom >*/
            int residue_sequence_number_;           /*!< Sequence number of the residue that the removed/replaced hydrogen atom belongs to >*/
            char residue_chain_id_;                 /*!< Chain id of the residue that the removed/replaced hydrogen atom belongs to >*/
            char residue_insertion_code_;           /*!< Insertion code of the residue that the removed/replaced hydrogen atom belongs to >*/
            char residue_alternate_location_;       /*!< Alternate location of the residue that the removed/replaced hydrogen atom belongs to >*/

    };
}


#endif // PDBPREPROCESSORREPLACEDHYDROGEN_HPP
