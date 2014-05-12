#ifndef PDBPREPROCESSORMISSINGRESIDUE_HPP
#define PDBPREPROCESSORMISSINGRESIDUE_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../../common.hpp"

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
            PdbPreprocessorMissingResidue(char chain_id, int start_index, int end_index, int index_before_gap, int index_after_gap, char starting_residue_insertion_code, char ending_residue_insertion_code);

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
              * An accessor function in order to access to the ending residue sequence number
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
              * An accessor function in order to access to the selected n termination
              * @return selected_n_termination_ attribute of the current object of this class
              */
            gmml::PossibleNChainTermination GetSelectedNTermination();
            /*! \fn
              * An accessor function in order to access to the selected c termination
              * @return selected_c_termination_ attribute of the current object of this class
              */
            gmml::PossibleCChainTermination GetSelectedCTermination();
            /*! \fn
              * An accessor function in order to access to the starting residue insertion code
              * @return starting_residue_insertion_code_ attribute of the current object of this class
              */
            char GetStartingResidueInsertionCode();
            /*! \fn
              * An accessor function in order to access to the ending residue insertion code
              * @return ending_residue_insertion_code_ attribute of the current object of this class
              */
            char GetEndingResidueInsertionCode();

            std::string GetStringFormatOfSelectedNTermination();
            std::string GetStringFormatOfSelectedCTermination();

            std::string GetStringFormatOfNTermination(gmml::PossibleNChainTermination n_termination);
            std::string GetStringFormatOfCTermination(gmml::PossibleCChainTermination c_termination);

            std::vector<std::string> GetAllPossibleNChainTerminationAsString();
            std::vector<std::string> GetAllPossibleCChainTerminationAsString();

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
              * A mutator function in order to set the selected n termination of the current object
              * Set the selected_n_termination_ attribute of the current pdb preprocessor missing residue
              * @param selected_n_termination The selected n termination attribute of the current object
              */
            void SetSelectedNTermination(gmml::PossibleNChainTermination selected_n_termination);
            /*! \fn
              * A mutator function in order to set the selected c termination of the current object
              * Set the selected_c_termination_ attribute of the current pdb preprocessor missing residue
              * @param selected_c_termination The selected c termination attribute of the current object
              */
            void SetSelectedCTermination(gmml::PossibleCChainTermination selected_c_termination);
            /*! \fn
              * A mutator function in order to set the starting residue insertion code of the current object
              * Set the starting_residue_insertion_code_ attribute of the current pdb preprocessor missing residue
              * @param starting_residue_insertion_code The starting residue insertion code attribute of the current object
              */
            void SetStartingResidueInsertionCode(char starting_residue_insertion_code);
            /*! \fn
              * A mutator function in order to set the ending residue insertion code of the current object
              * Set the ending_residue_insertion_code_ attribute of the current pdb preprocessor missing residue
              * @param ending_residue_insertion_code The ending residue insertion code attribute of the current object
              */
            void SetEndingResidueInsertionCode(char ending_residue_insertion_code);

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
            gmml::PossibleNChainTermination selected_n_termination_;
            gmml::PossibleCChainTermination selected_c_termination_;
            char starting_residue_insertion_code_;
            char ending_residue_insertion_code_;

    };
}


#endif // PDBPREPROCESSORMISSINGRESIDUE_HPP
