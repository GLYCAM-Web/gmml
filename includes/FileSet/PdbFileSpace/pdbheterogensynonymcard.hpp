// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGENSYNONYMCARD_HPP
#define PDBHETEROGENSYNONYMCARD_HPP

#include <string>
#include <map>
#include <sstream>

namespace PdbFileSpace
{
    class PdbHeterogenSynonym;
    class PdbHeterogenSynonymCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, PdbHeterogenSynonym*> HeterogenSynonymMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHeterogenSynonymCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbHeterogenSynonymCard(const std::string& record_name);
            PdbHeterogenSynonymCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a heterogen synonym card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the heterogen synonym in a heterogen synonym card
              * @return heterogen_synonym_ attribute of the current object of this class
              */
            HeterogenSynonymMap GetHeterogensSynonyms();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current heterogen synonym card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);

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
            std::string record_name_;
            HeterogenSynonymMap heterogens_synonyms_;
    };
}
#endif // PDBHETEROGENSYNONYMCARD_HPP
