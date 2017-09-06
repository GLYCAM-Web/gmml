// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia// Author: Alireza Khatamian

#ifndef PDBRESIDUEMODIFICATIONCARD_HPP
#define PDBRESIDUEMODIFICATIONCARD_HPP

#include <map>
#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbResidueModification;
    class PdbResidueModificationCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between residue modification id code and the residue modification itself
              */
            typedef std::map<std::string, PdbResidueModification*> ResidueModificationMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbResidueModificationCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbResidueModificationCard(const std::string& record_name);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbResidueModificationCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a residue modification card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the residue modifications in a residue modification card
              * @return residue_modifications_ attribute of the current object of this class
              */
            ResidueModificationMap GetResidueModifications();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current residue modification card
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
              * A function to print out the pdb residue modification card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;                           /*!< Record name of this card appears in the first column of each line in a pdb file >*/
            ResidueModificationMap residue_modifications_;      /*!< Residue modification object map consists of all residue modification objects in the residue modification card of a pdb file >*/
    };
}

#endif // PDBRESIDUEMODIFICATIONCARD_HPP
