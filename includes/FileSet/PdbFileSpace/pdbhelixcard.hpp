// Author: Alireza Khatamian

#ifndef PDBHELIXCARD_HPP
#define PDBHELIXCARD_HPP

#include <string>
#include <map>

namespace PdbFileSpace
{
    class PdbHelix;
    class PdbHelixCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, PdbHelix*> HelixMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbHelixCard();
            PdbHelixCard(const std::string& record_name);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            HelixMap GetHelixes();

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
            HelixMap helixes_;
    };
}
#endif // PDBHELIXCARD_HPP
