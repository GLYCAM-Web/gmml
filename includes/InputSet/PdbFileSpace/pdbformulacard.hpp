// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBFORMULACARD_HPP
#define PDBFORMULACARD_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbFormulaCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbFormulaCard();
            /*! \fn
              * Constructor with required parameters
              * @param heterogen_identifier Heterogen identifier
              * @param component_number Component number of the formula
              * @param chemical_formula Chemical formula of the object
              */
            PdbFormulaCard(const std::string& heterogen_identifier, int component_number, const std::string& chemical_formula);
            /*! \fn
              * Constructor with required parameters
              * @param specification_block
              */
            PdbFormulaCard(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */

            /*! \fn
              * An accessor function in order to access to the heterogen identifier in a pdb formula
              * @return heterogen_identifier_ attribute of the current object of this class
              */
            std::string GetHeterogenIdentifier();
            /*! \fn
              * An accessor function in order to access to the component number in a pdb formula
              * @return component_number_ attribute of the current object of this class
              */
            int GetComponentNumber();
            /*! \fn
              * An accessor function in order to access to the chemical formula in a pdb formula
              * @return chemical_formula_ attribute of the current object of this class
              */
            std::string GetChemicalFormula();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //

           //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the heterogen identifier of the current object
              * Set the heterogen_identifier_ attribute of the current formula
              * @param heterogen_identifier The heterogen identifier of the current object
              */
            void SetHeterogenIdentifier(const std::string heterogen_identifier);
            /*! \fn
              * A mutator function in order to set the component number of the current object
              * Set the component_number_ attribute of the current formula
              * @param component_number The component number of the current object
              */
            void SetComponentNumber(int component_number);
            /*! \fn
              * A mutator function in order to set the chemical formula of the current object
              * Set the chemical_formula_ attribute of the current formula
              * @param chemical_formula The chemical formula of the current object
              */
            void SetChemicalFormula(const std::string chemical_formula);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the formula contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string heterogen_identifier_;      /*!< Heterogen identifier >*/
            int component_number_;                  /*!< Component number >*/
            std::string chemical_formula_;          /*!< Chemical formula of the object >*/
    };
}
#endif // PDBFORMULACARD_HPP
