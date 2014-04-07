#ifndef PDBPREPROCESSORCHAINTERMINATION_HPP
#define PDBPREPROCESSORCHAINTERMINATION_HPP

#include <string>
#include <iostream>
#include <vector>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorChainTermination
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorChainTermination();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue chain Id
              * @return residue_chain_id_ attribute of the current object of this class
              */
            char GetResidueChainId();
            /*! \fn
              * An accessor function in order to access to the residue start index
              * @return residue_start_index_ attribute of the current object of this class
              */
            int GetResidueStartIndex();
            /*! \fn
              * An accessor function in order to access to the residue end index
              * @return residue_end_index_ attribute of the current object of this class
              */
            int GetResidueEndIndex();
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
              * An accessor function in order to access to the selected n terminations
              * @return selected_n_terminations_ attribute of the current object of this class
              */
            std::string GetSelectedNTerminations();
            /*! \fn
              * An accessor function in order to access to the selected c terminations
              * @return selected_c_terminations_ attribute of the current object of this class
              */
            std::string GetSelectedCTerminations();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current pdb preprocessor chain termination
              * @param residue_chain_id The residue chain id attribute of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the residue start index of the current object
              * Set the residue_start_index_ attribute of the current pdb preprocessor chain termination
              * @param residue_start_index The residue start index attribute of the current object
              */
            void SetResidueStartIndex(int residue_start_index);
            /*! \fn
              * A mutator function in order to set the residue end index of the current object
              * Set the residue_end_index_ attribute of the current pdb preprocessor chain termination
              * @param residue_end_index The residue end index attribute of the current object
              */
            void SetResidueEndIndex(int residue_end_index);
            /*! \fn
              * A mutator function in order to set the possible n terminations of the current object
              * Set the possible_n_terminations_ attribute of the current pdb preprocessor chain termination
              * @param possible_n_terminations The possible n terminations attribute of the current object
              */
            void SetPossibleNTerminations(std::vector<std::string> possible_n_terminations);
            /*! \fn
              * A mutator function in order to set the possible c terminations of the current object
              * Set the possible_c_terminations_ attribute of the current pdb preprocessor chain termination
              * @param possible_c_terminations The possible c terminations attribute of the current object
              */
            void SetPossibleCTerminations(std::vector<std::string> possible_c_terminations);
            /*! \fn
              * A mutator function in order to set the selected n terminations of the current object
              * Set the selected_n_terminations_ attribute of the current pdb preprocessor chain termination
              * @param selected_n_terminations The selected n terminations attribute of the current object
              */
            void SetSelectedNTerminations(std::string selected_n_terminations);
            /*! \fn
              * A mutator function in order to set the selected c terminations of the current object
              * Set the selected_c_terminations_ attribute of the current pdb preprocessor chain termination
              * @param selected_c_terminations The selected c terminations attribute of the current object
              */
            void SetSelectedCTerminations(std::string selected_c_terminations);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor chain termination contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_;
            int residue_start_index_;
            int residue_end_index_;
            std::vector<std::string> possible_n_terminations_;
            std::vector<std::string> possible_c_terminations_;
            std::string selected_n_terminations_;
            std::string selected_c_terminations_;

    };
}

#endif // PDBPREPROCESSORCHAINTERMINATION_HPP
