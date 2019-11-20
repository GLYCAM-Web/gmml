// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBHETEROGENNAMESECTION_HPP
#define PDBHETEROGENNAMESECTION_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbHeterogenNameCard;
    class PdbHeterogenNameSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between heterogen identifier and heterogen name
              */
            typedef std::map<std::string, PdbHeterogenNameCard*> HeterogenNameCardMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHeterogenNameSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbHeterogenNameSection(const std::string& record_name);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbHeterogenNameSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a heterogen name card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the heterogen names in a heterogen name card
              * @return heterogen_name_cards_ attribute of the current object of this class
              */
            HeterogenNameCardMap GetHeterogenNameCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current heterogen name card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
/** @*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the heterogen name card contents in a structural format
              * Print out the information in a deheterogen_cards_name_fined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of heterogen name card in a pdb file >*/
            HeterogenNameCardMap heterogen_name_cards_;  /*!< Heterogen name map >*/
    };
}

#endif // PDBHETEROGENNAMESECTION_HPP
