// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBRESIDUESEQUENCECARD_HPP
#define PDBRESIDUESEQUENCECARD_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbResidueSequenceCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbResidueSequenceCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name Record name which is the first column of each line in a pdb file
              * @param number_of_residues  Number of residues involving in a sequence residue
              * @param residue_names Name of residues in the sequence residue
              */
            PdbResidueSequenceCard(char chain_id, int number_of_residues, const std::vector<std::string>& residue_names);
            /*! \fn
              * Constructor with required parameters
              * @param specification_block
              */
            PdbResidueSequenceCard(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the chain id in residue sequence
              * @return chain_id_ attribute of the current object of this class
              */
            char GetChainId();
            /*! \fn
              * An accessor function in order to access to the number of residues in residue sequence
              * @return number_of_residues_ attribute of the current object of this class
              */
            int GetNumberOfResidues();
            /*! \fn
              * An accessor function in order to access to the list of residue names in residue sequence
              * @return residue_names_ of the current object of this class
              */
            std::vector<std::string> GetResidueNames();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the chain id of the current object
              * Set the chain_id_ attribute of the current residue sequence
              * @param chain_id The chain id of the current object
              */
            void SetChainId(char chain_id);
            /*! \fn
              * A mutator function in order to set the number of residues of the current object
              * Set number_of_residues_ attribute of the current residue sequence
              * @param number_of_residues The number of residues of the current object
              */
            void SetNumberOfResidues(int number_of_residues);
            /*! \fn
              * A mutator function in order to set the list of residue names of the current object
              * Set the residue_names_ of the current residue sequence
              * @param residue_names The residue names of the current object
              */
            void SetResidueNames(const std::vector<std::string> residue_names);
            /*! \fn
              * A function in order to add the residue name to the current object
              * Set the residue_name_ attribute of the current residue sequence
              * @param residue_name The residue name of the current object
              */
            void AddResidueName(const std::string residue_name);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb residue sequence contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            char chain_id_;                         /*!< Chain identifier of the residue sequence >*/
            int number_of_residues_;                /*!< Number of residues in the residue sequence >*/
            std::vector<std::string> residue_names_;/*!< Name of all residues in the residue sequence >*/
    };
}

#endif // PDBRESIDUESEQUENCECARD_HPP
