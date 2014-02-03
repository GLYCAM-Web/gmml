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
            PdbFormula();
            PdbFormula(const std::string& heterogen_identifier, int component_number, const std::string& chemical_formula);
            PdbFormula(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetHeterogenIdentifier();
            int GetComponentNumber();
            std::string GetChemicalFormula();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetHeterogenIdentifier(const std::string heterogen_identifier);
            void SetComponentNumber(int component_number);
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
