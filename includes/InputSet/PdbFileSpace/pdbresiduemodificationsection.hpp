// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBRESIDUEMODIFICATIONSECTION_HPP
#define PDBRESIDUEMODIFICATIONSECTION_HPP

#include <map>
#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbResidueModificationCard;
    class PdbResidueModificationSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between residue modification id code and the residue modification itself
              */
            typedef std::map<std::string, PdbResidueModificationCard*> ResidueModificationCardMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbResidueModificationSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbResidueModificationSection(const std::string& record_name);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbResidueModificationSection(std::stringstream& stream_block);

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
              * @return residue_modification_cards_ attribute of the current object of this class
              */
            ResidueModificationCardMap GetResidueModificationCards();
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
              * @param record_name The record name PdbResidueModificationSectionof the current object
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
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;                           /*!< Record name of this card appears in the first column of each line in a pdb file >*/
            ResidueModificationCardMap residue_modification_cards_;      /*!< Residue modification object map consists of all residue modification objects in the residue modification card of a pdb file >*/
    };
}

#endif // PDBRESIDUEMODIFICATIONSECTION_HPP
