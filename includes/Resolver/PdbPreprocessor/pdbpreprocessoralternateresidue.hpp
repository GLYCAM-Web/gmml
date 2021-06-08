#ifndef PDBPREPROCESSORALTERNATERESIDUE_HPP
#define PDBPREPROCESSORALTERNATERESIDUE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorAlternateResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorAlternateResidue();

            /*! \fn
              * Constructor to initialize the attributes of an entity of an alternate resiude
              * (Residues with same description and different alternate location attributes are declared as alternate residues)
              * @param residue_name Name of residue that has alternate location(s)
              * @param chain_id Chain id of the residue with alternate location(s)
              * @param sequence_number Sequence number of the residue with alternate location(s)
              * @param residue_insertion_code Insertion code of the residue with alternate location(s)
              * @param residue_alternate_location List of alternate locations that the residue has
              * @param selected_alternate_location List of true/false values corresponding to each alternate location that indicates which one is selected
              */
            PdbPreprocessorAlternateResidue(std::string residue_name, char chain_id, int sequence_number, char residue_insertion_code,
                                            std::vector<char> residue_alternate_location, std::vector<bool> selected_alternate_location);

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
            std::vector<char> GetResidueAlternateLocation();
            /*! \fn
              * An accessor function in order to access to the selected alternate location
              * @return selected_alternate_location_ attribute of the current object of this class
              */
            std::vector<bool> GetSelectedAlternateLocation();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current pdb preprocessor alternate residue
              * @param residue_chain_id The residue chain id attribute of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the residue sequence number of the current object
              * Set the residue_sequence_number_ attribute of the current pdb preprocessor alternate residue
              * @param residue_sequence_number The residue sequence number attribute of the current object
              */
            void SetResidueSequenceNumber(int residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current pdb preprocessor alternate residue
              * @param residue_name The residue name attribute of the current object
              */
            void SetResidueName(std::string residue_name);
            /*! \fn
              * A mutator function in order to set the residue insertion code of the current object
              * Set the residue_insertion_code_ attribute of the current pdb preprocessor alternate residue
              * @param residue_insertion_code The residue insertion code attribute of the current object
              */
            void SetResidueInsertionCode(char residue_insertion_code);
            /*! \fn
              * A mutator function in order to set the residue alternate location of the current object
              * Set the residue_alternate_location_ attribute of the current pdb preprocessor alternate residue
              * @param residue_alternate_location The residue alternate location attribute of the current object
              */
            void SetResidueAlternateLocation(std::vector<char> residue_alternate_location);
            /*! \fn
              * A mutator function in order to set the selected alternate location of the current object
              * Set the selected_alternate_location_ attribute of the current pdb preprocessor alternate residue
              * @param selected_alternate_location The selected alternate location attribute of the current object
              */
            void SetSelectedAlternateLocation(std::vector<bool> selected_alternate_location);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor alternate residue in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_;                             /*!< Chain id of a residue with alternate location(s) >*/
            int residue_sequence_number_;                       /*!< Sequence number of a residue with alternate location(s) >*/
            std::string residue_name_;                          /*!< Name of the residue with alternate location(s) >*/
            char residue_insertion_code_;                       /*!< Insertion code of the residue with alternate location(s) >*/
            std::vector<char> residue_alternate_location_;      /*!< List of alternate locations of the residue with different possible location(s) >*/
            std::vector<bool> selected_alternate_location_;     /*!< List of true/false values that indicate which one of the alternate location(s) is selected >*/

    };
}

#endif // PDBPREPROCESSORALTERNATERESIDUE_HPP
