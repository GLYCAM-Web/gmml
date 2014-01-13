// Author: Alireza Khatamian

#ifndef PDBHETEROGENCARD_HPP
#define PDBHETEROGENCARD_HPP

#include <string>
#include <map>
#include <sstream>

namespace PdbFileSpace
{
    class PdbHeterogen;
    class PdbHeterogenCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, PdbHeterogen*> HeterogenMap;
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbHeterogenCard();
            PdbHeterogenCard(const std::string& record_name);
            PdbHeterogenCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            HeterogenMap GetHeterogens();

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
            HeterogenMap heterogens_;
    };
}

#endif // PDBHETEROGENCARD_HPP
