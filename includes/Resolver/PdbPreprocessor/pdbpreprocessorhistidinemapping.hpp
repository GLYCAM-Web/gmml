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
              * An accessor function in order to access to the residue number
              * @return residue_number_ attribute of the current object of this class
              */
            int GetResidueNumber();
            /*! \fn
              * An accessor function in order to access to the possible_mappings
              * @return possible_mappings_ attribute of the current object of this class
              */
            std::vector<std::string> GetPossibleMappings();
            /*! \fn
              * An accessor function in order to access to the selected mappings
              * @return selected_mappings_ attribute of the current object of this class
              */
            std::string GetSelectedMappings();

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
              * A mutator function in order to set the residue number of the current object
              * Set the residue_number_ attribute of the current pdb preprocessor histidine mapping
              * @param residue_number The residue number attribute of the current object
              */
            void SetResidueNumber(int residue_number);
            /*! \fn
              * A mutator function in order to set the possible_mappings of the current object
              * Set the possible_mappings_ attribute of the current pdb preprocessor histidine mapping
              * @param possible_mappings The possible mappings attribute of the current object
              */
            void SetPossibleMappings(std::vector<std::string> possible_mappings);
            /*! \fn
              * A mutator function in order to set the selected mappings of the current object
              * Set the selected_mappings_ attribute of the current pdb preprocessor histidine mapping
              * @param selected_mappings The selected mappings attribute of the current object
              */
            void SetSelectedMappings(std::string selected_mappings);

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
            int residue_number_;
            std::vector<std::string> possible_mappings_;
            std::string selected_mappings_;

    };
}


#endif // PDBPREPROCESSORHISTIDINEMAPPING_HPP
