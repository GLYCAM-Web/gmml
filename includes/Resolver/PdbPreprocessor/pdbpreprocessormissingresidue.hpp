#ifndef PDBPREPROCESSORMISSINGRESIDUE_HPP
#define PDBPREPROCESSORMISSINGRESIDUE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorMissingResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorMissingResidue();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue chain Id
              * @return residue_chain_id_ attribute of the current object of this class
              */
            char GetResidueChainId();
            /*! \fn
              * An accessor function in order to access to the starting residue sequence number
              * @return starting_residue_sequence_number_ attribute of the current object of this class
              */
            int GetStartingResidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the starting residue sequence number
              * @return ending_residue_sequence_number_ attribute of the current object of this class
              */
            int GetEndingResidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the residue before gap
              * @return residue_before_gap_ attribute of the current object of this class
              */
            int GetResidueBeforeGap();
            /*! \fn
              * An accessor function in order to access to the residue after gap
              * @return residue_after_gap_ attribute of the current object of this class
              */
            int GetResidueAfterGap();
            /*! \fn
              * An accessor function in order to access to the possible n teminations
              * @return possible_n_terminations_ attribute of the current object of this class
              */
            std::vector<std::string> GetPossibleNTerminations();
            /*! \fn
              * An accessor function in order to access to the possible c teminations
              * @return possible_c_terminations_ attribute of the current object of this class
              */
            std::vector<std::string> GetPossibleCTerminations();
            /*! \fn
              * An accessor function in order to access to the selected n termination
              * @return selected_n_termination_ attribute of the current object of this class
              */
            std::string GetSelectedNTermination();
            /*! \fn
              * An accessor function in order to access to the selected c termination
              * @return selected_c_termination_ attribute of the current object of this class
              */
            std::string GetSelectedCTermination();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current pdb preprocessor missing residue
              * @param residue_chain_id The residue chain id attribute of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the starting residue sequence number of the current object
              * Set the starting_residue_sequence_number_ attribute of the current pdb preprocessor missing residue
              * @param starting_residue_sequence_number The starting residue sequence number attribute of the current object
              */
            void SetStartingResidueSequenceNumber(int starting_residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the ending residue sequence number of the current object
              * Set the ending_residue_sequence_number_ attribute of the current pdb preprocessor missing residue
              * @param ending_residue_sequence_number The ending residue sequence number attribute of the current object
              */
            void SetEndingResidueSequenceNumber(int ending_residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the residue_before_gap of the current object
              * Set the residue_before_gap_ attribute of the current pdb preprocessor missing residue
              * @param residue_before_gap The residue before gap attribute of the current object
              */
            void SetResidueBeforeGap(int residue_before_gap);
            /*! \fn
              * A mutator function in order to set the residue after gap of the current object
              * Set the residue_after_gap_ attribute of the current pdb preprocessor missing residue
              * @param residue_after_gap The residue after gap attribute of the current object
              */
            void SetResidueAfterGap(int residue_after_gap);
            /*! \fn
              * A mutator function in order to set the possible n terminations of the current object
              * Set the possible_n_terminations_ attribute of the current pdb preprocessor missing residue
              * @param possible_n_terminations The possible n terminations attribute of the current object
              */
            void SetPossibleNTerminations(std::vector<std::string> possible_n_terminations);
            /*! \fn
              * A mutator function in order to set the possible c terminations of the current object
              * Set the possible_c_terminations_ attribute of the current pdb preprocessor missing residue
              * @param possible_c_terminations The possible c terminations attribute of the current object
              */
            void SetPossibleCTerminations(std::vector<std::string> possible_c_terminations);
            /*! \fn
              * A mutator function in order to set the selected n termination of the current object
              * Set the selected_n_termination_ attribute of the current pdb preprocessor missing residue
              * @param selected_n_termination The selected n termination attribute of the current object
              */
            void SetSelectedNTermination(std::string selected_n_termination);
            /*! \fn
              * A mutator function in order to set the selected c termination of the current object
              * Set the selected_c_termination_ attribute of the current pdb preprocessor missing residue
              * @param selected_c_termination The selected c termination attribute of the current object
              */
            void SetSelectedCTermination(std::string selected_c_termination);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor missing residue contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_;
            int starting_residue_sequence_number_;
            int ending_residue_sequence_number_;
            int residue_before_gap_;
            int residue_after_gap_;
            std::vector<std::string> possible_n_terminations_;
            std::vector<std::string> possible_c_terminations_;
            std::string selected_n_termination_;
            std::string selected_c_termination_;

    };
}


#endif // PDBPREPROCESSORMISSINGRESIDUE_HPP
