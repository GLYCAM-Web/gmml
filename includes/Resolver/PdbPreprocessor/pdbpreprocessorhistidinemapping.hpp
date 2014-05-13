#ifndef PDBPREPROCESSORHISTIDINEMAPPING_HPP
#define PDBPREPROCESSORHISTIDINEMAPPING_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../../common.hpp"

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

            PdbPreprocessorHistidineMapping(char chain_id, int residue_sequence_number, gmml::PdbPreprocessorHISMapping selected_mapping, char residue_insertion_code);

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
              * An accessor function in order to access to the selected mapping
              * @return selected_mapping_ attribute of the current object of this class
              */
            gmml::PdbPreprocessorHISMapping GetSelectedMapping();
            /*! \fn
              * An accessor function in order to access to the residue insertion code
              * @return residue_insertion_code_ attribute of the current object of this class
              */
            char GetResidueInsertionCode();

            /*! \fn
              * An accessor function in order to access to the string format of selected mapping
              * @return
              */
            std::string GetStringFormatOfSelectedMapping();
            /*! \fn
              * An accessor function in order to access to the string format of a selected mapping
              * @param mapping The histidine mapping of the current object
              * @return
              */
            std::string GetStringFormatOfMapping(gmml::PdbPreprocessorHISMapping mapping);
            /*! \fn
              * An accessor function in order to access to all HIS mappings as string
              * @return all_his_mapping_as_string
              */
            std::vector<std::string> GetAllHISMappingAsString();

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
              * A mutator function in order to set the selected mapping of the current object
              * Set the selected_mapping_ attribute of the current pdb preprocessor histidine mapping
              * @param selected_mapping The selected mapping attribute of the current object
              */
            void SetSelectedMapping(gmml::PdbPreprocessorHISMapping selected_mapping);
            /*! \fn
              * A mutator function in order to set the insertion code of the current object
              * Set the insertion_code_ attribute of the current pdb preprocessor histidine mapping
              * @param insertion_code The insertion code attribute of the current object
              */
            void SetResidueInsertionCode(char residue_insertion_code);

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
            gmml::PdbPreprocessorHISMapping selected_mapping_;
            char residue_insertion_code_;

    };
}


#endif // PDBPREPROCESSORHISTIDINEMAPPING_HPP
