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

            /*! \fn
              * Constructor to initialized the attributes of an entity of histidine residue
              * Residues with HIS name are defined as histidine residue
              * @param chain_id Chain id of the histidine residue
              * @param residue_sequence_number Sequence number of the histidine residue
              * @param selected_mapping Selected type of the histidine residue from the all possible types that have been defined in PdbPreprocessorHISMapping enumerator
              * @param residue_insertion_code Insertion code of the histidine residue
              * @param residue_alternate_location Alternate location of the histidine residue
              */
            PdbPreprocessorHistidineMapping(char chain_id, int residue_sequence_number, gmml::PdbPreprocessorHISMapping selected_mapping,
                                            char residue_insertion_code, char residue_alternate_location);

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
              * An accessor function in order to access to the residue alternate location
              * @return residue_alternate_location_ attribute of the current object of this class
              */
            char GetResidueAlternateLocation();

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
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
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
              * Set the residue_insertion_code_ attribute of the current pdb preprocessor histidine mapping
              * @param insertion_code The residue insertion code attribute of the current object
              */
            void SetResidueInsertionCode(char residue_insertion_code);
            /*! \fn
              * A mutator function in order to set the residue alternate location of the current object
              * Set the residue_alternate_location_ attribute of the current pdb preprocessor histidine mapping
              * @param residue_alternate_location The residue alternate location attribute of the current object
              */
            void SetResidueAlternateLocation(char residue_alternate_location);
/** @*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor histidine mapping contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_;                                 /*!< Chain id of the histidine resiude >*/
            int residue_sequence_number_;                           /*!< Sequence number of the histidine residue >*/
            gmml::PdbPreprocessorHISMapping selected_mapping_;      /*!< Selected mapping of the histidine residue >*/
            char residue_insertion_code_;                           /*!< Insertion code of the histidine residue >*/
            char residue_alternate_location_;                       /*!< Alternate location of the histidine residue >*/

    };
}


#endif // PDBPREPROCESSORHISTIDINEMAPPING_HPP
