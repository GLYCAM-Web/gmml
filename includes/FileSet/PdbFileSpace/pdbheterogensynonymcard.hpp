// Author: Alireza Khatamian

#ifndef PDBHETEROGENSYNONYMCARD_HPP
#define PDBHETEROGENSYNONYMCARD_HPP

#include <string>
#include <map>

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
            PdbHeterogenSynonymCard();
            PdbHeterogenSynonymCard(const std::string& record_name);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            HeterogenSynonymMap GetHeterogensSynonyms();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
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
