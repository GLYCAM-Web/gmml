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
            /*! \fn
              * Default constructor
              */
            PdbHeterogenSynonym();
            /*! \fn
              * Constructor with required parameters
              * @param heterogen_identifier
              * @param heterogen_synonyms
              */
            PdbHeterogenSynonym(const std::string& heterogen_identifier, const std::vector<std::string>& heterogen_synonyms);
            PdbHeterogenSynonym(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the heterogen identifier in a heterogen synonym
              * @return heterogen_identifier_ attribute of the current object of this class
              */
            std::string GetHeterogenIdentifier();
            /*! \fn
              * An accessor function in order to access to the list of heterogen synonyms in a heterogen synonym
              * @return heterogen_synonyms_ of the current object of this class
              */
            std::vector<std::string> GetHeterogenSynonyms();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the heterogen identifier of the current object
              * Set the heterogen_identifier_ attribute of the current heterogen synonym
              * @param heterogen_identifier The heterogen identifier of the current object
              */
            void SetHeterogenIdentifier(const std::string heterogen_identifier);
            /*! \fn
              * A mutator function in order to set the list of heterogen synonyms of the current object
              * Set the heterogen_synonyms_ of the current heterogen synonym
              * @param heterogen_synonyms The heterogen synonyms of the current object
              */
            void SetHeterogenSynonyms(const std::vector<std::string> heterogen_synonyms);
            /*! \fn
              * A function in order to add the heterogen synonyms to the current object
              * Set the heterogen_synonyms_ of the current heterogen synonym
              * @param heterogen_synonyms The heterogen synonyms of the current object
              */
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
