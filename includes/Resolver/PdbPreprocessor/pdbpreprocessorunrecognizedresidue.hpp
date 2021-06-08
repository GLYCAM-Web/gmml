#ifndef PDBPREPROCESSORUNRECOGNIZEDRESIDUE_HPP
#define PDBPREPROCESSORUNRECOGNIZEDRESIDUE_HPP

#include <string>
#include <iostream>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorUnrecognizedResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorUnrecognizedResidue();

            /*! \fn
              * Constructor to initialized the attributes of the object
              * @param residue_name Name of the unrecognized residue
              * @param chain_id Chain id of the unrecognized residue
              * @param sequence_number Sequence number of the unrecognnized residue
              * @param residue_insertion_code Insertion code of the unrecognized residue
              * @param residue_alternate_location Alternated location of the unrecognized residue
              */
            PdbPreprocessorUnrecognizedResidue(
                std::string residue_name, char chain_id, int sequence_number,
                char residue_insertion_code, char residue_alternate_location,
                bool middle_of_chain);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the residue chain Id
              * @return residue_chain_id_ attribute of the current object of this class
              */
            char GetResidueChainId();
            /*! \fn
              * An accessor function in order to access to the residue sequence number
              * @return residue_sequence_number_ attribute of the current object of this class
              */
            int GetResidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the residue name
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
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
            /*! \fn
              * An accessor function in order to access to the middle of chain attribute of the current object
              * @return middle_of_chain_ attribute of the current object of this class
              */
            bool GetMiddleOfChain();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current pdb preprocessor unrecognized residue
              * @param residue_chain_id The residue chain id attribute of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the residue sequence number of the current object
              * Set the residue_sequence_number_ attribute of the current pdb preprocessor unrecognized residue
              * @param residue_sequence_number The residue sequence number attribute of the current object
              */
            void SetResidueSequenceNumber(int residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current pdb preprocessor unrecognized residue
              * @param residue_name The residue name attribute of the current object
              */
            void SetResidueName(std::string residue_name);
            /*! \fn
              * A mutator function in order to set the residue insertion code of the current object
              * Set the residue_insertion_code_ attribute of the current pdb preprocessor unrecognized residue
              * @param residue_insertion_code The residue insertion code attribute of the current object
              */
            void SetResidueInsertionCode(char residue_insertion_code);
            /*! \fn
              * A mutator function in order to set the residue alternate location of the current object
              * Set the residue_alternate_location_ attribute of the current pdb preprocessor unrecognized residue
              * @param residue_alternate_location The residue alternate location attribute of the current object
              */
            void SetResidueAlternateLocation(char residue_alternate_location);
            /*! \fn
              * A mutator function in order to set if this residue is in the middle of a chain in the current object
              * Set the middle_of_chain_ attribute of the current pdb preprocessor unrecognized residue
              * @param middle_of_chain The boolean value of middle of chain attribute of the current object
              */
            void SetMiddleOfChain(bool middle_of_chain);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor unrecognized residue contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_;                     /*!< Chain id of the unrecognized residue in a pdb file >*/
            int residue_sequence_number_;               /*!< Sequence number of the unrecognized residue in an specific chain id >*/
            std::string residue_name_;                  /*!< Name of the unrecognized residue >*/
            char residue_insertion_code_;               /*!< Insertion code of the unrecognized residue >*/
            char residue_alternate_location_;           /*!< Alternate location of the unrecognized residue >*/
            bool middle_of_chain_;
    };
}

#endif // PDBPREPROCESSORUNRECOGNIZEDRESIDUE_HPP
