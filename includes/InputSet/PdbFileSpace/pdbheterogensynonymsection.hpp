// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBHETEROGENSYNONYMSECTION_HPP
#define PDBHETEROGENSYNONYMSECTION_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbHeterogenSynonymCard;
    class PdbHeterogenSynonymSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between heterogen identifier and heterogen synonym
              */
            typedef std::map<std::string, PdbHeterogenSynonymCard*> HeterogenSynonymCardMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHeterogenSynonymSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbHeterogenSynonymSection(const std::string& record_name);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbHeterogenSynonymSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a heterogen synonym card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the heterogen synonym in a heterogen synonym card
              * @return heterogen_synonym_ attribute of the current object of this class
              */
            HeterogenSynonymCardMap GetHeterogensSynonymCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current heterogen synonym card
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
              * A function to print out the heterogen synonym card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;                   /*!< Record name of a heterogen synonym card >*/
            HeterogenSynonymCardMap heterogens_synonym_cards_;   /*!< Map of heterogen synonyms with heterogen identifier as key >*/
    };
}
#endif // PDBHETEROGENSYNONYMSECTION_HPP
