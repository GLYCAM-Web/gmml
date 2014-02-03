// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBFORMULA_HPP
#define PDBFORMULA_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbFormula
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
            PdbFormula();
            /*! \fn
              * Constructor with required parameters
              * @param heterogen_identifier
              * @param component_number
              * @param chemical_formula
              */
            PdbFormula(const std::string& heterogen_identifier, int component_number, const std::string& chemical_formula);
            PdbFormula(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string heterogen_identifier_;
            int component_number_;
            std::string chemical_formula_;
    };
}
#endif // PDBFORMULA_HPP
