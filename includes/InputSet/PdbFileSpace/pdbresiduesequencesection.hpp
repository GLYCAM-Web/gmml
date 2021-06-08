// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBRESIDUESEQUENCESECTION_HPP
#define PDBRESIDUESEQUENCESECTION_HPP

#include <string>
#include <map>
#include <iostream>

namespace PdbFileSpace
{
    class PdbResidueSequenceCard;
    class PdbResidueSequenceSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between chain identifier of a residue sequence and the residue sequence itself
              */
            typedef std::map<char, PdbResidueSequenceCard*> ResidueSequenceCardMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbResidueSequenceSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbResidueSequenceSection(const std::string& record_name);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbResidueSequenceSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in residue sequence card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the residue sequence chains in residue sequence card
              * @return residue_sequence_chains_ attribute of the current object of this class
              */
            ResidueSequenceCardMap GetResidueSequenceChain();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current residue sequence card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb residue sequence card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;                       /*!< Name of the residue sequence record which is the first column in each line of a pdb file >*/
            ResidueSequenceCardMap residue_sequence_chains_;    /*!< Mapping of all residue sequences with chain identifier as the key >*/

    };
}

#endif // PDBRESIDUESEQUENCESECTION_HPP
