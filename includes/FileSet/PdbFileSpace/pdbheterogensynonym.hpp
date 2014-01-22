// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGENSYNONYM_HPP
#define PDBHETEROGENSYNONYM_HPP

#include <string>
#include <vector>

namespace PdbFileSpace
{
    class PdbHeterogenSynonym
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbHeterogenSynonym();
            PdbHeterogenSynonym(const std::string& heterogen_identifier, const std::vector<std::string>& heterogen_synonyms);
            PdbHeterogenSynonym(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetHeterogenIdentifier();
            std::vector<std::string> GetHeterogenSynonyms();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetHeterogenIdentifier(const std::string heterogen_identifier);
            void SetHeterogenSynonyms(const std::vector<std::string> heterogen_synonyms);
            void AddHeterogenSynonym(const std::string heterogen_synonym);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string heterogen_identifier_;
            std::vector<std::string> heterogen_synonyms_;
    };
}
#endif // PDBHETEROGENSYNONYM_HPP
