#ifndef PDBPREPROCESSORHISTIDINEMAPPING_HPP
#define PDBPREPROCESSORHISTIDINEMAPPING_HPP

#include <string>
#include <iostream>
#include <vector>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorHistidineMapping
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorHistidineMapping();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
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
              * An accessor function in order to access to the possible_mappings
              * @return possible_mappings_ attribute of the current object of this class
              */
            std::vector<std::string> GetPossibleMappings();
            /*! \fn
              * An accessor function in order to access to the selected mapping
              * @return selected_mapping_ attribute of the current object of this class
              */
            std::string GetSelectedMapping();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current pdb preprocessor histidine mapping
              * @param residue_chain_id The residue chain id attribute of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the residue sequence number of the current object
              * Set the residue_sequence_number_ attribute of the current pdb preprocessor histidine mapping
              * @param residue_sequence_number The residue sequence number attribute of the current object
              */
            void SetResidueSequenceNumber(int residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the possible_mappings of the current object
              * Set the possible_mappings_ attribute of the current pdb preprocessor histidine mapping
              * @param possible_mappings The possible mappings attribute of the current object
              */
            void SetPossibleMappings(std::vector<std::string> possible_mappings);
            /*! \fn
              * A mutator function in order to set the selected mapping of the current object
              * Set the selected_mapping_ attribute of the current pdb preprocessor histidine mapping
              * @param selected_mapping The selected mapping attribute of the current object
              */
            void SetSelectedMapping(std::string selected_mapping);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor histidine mapping contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_;
            int residue_sequence_number_;
            std::vector<std::string> possible_mappings_;
            std::string selected_mapping_;

    };
}


#endif // PDBPREPROCESSORHISTIDINEMAPPING_HPP
