// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBFORMULACARD_HPP
#define PDBFORMULACARD_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbFormula;
    class PdbFormulaCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between heterogen identifier and its chemical formula
              */
            typedef std::map<std::string, PdbFormula*> FormulaMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbFormulaCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbFormulaCard(const std::string& record_name);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbFormulaCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a formula card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the formulas in a formula card
              * @return formulas_ attribute of the current object of this class
              */
            FormulaMap GetFormulas();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current formula card
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
              * A function to print out the formula card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;       /*!< Record name of formula card in a pdb file >*/
            FormulaMap formulas_;           /*!< Map of formulas with heterogen identifier as key >*/
    };
}
#endif // PDBFORMULACARD_HPP
