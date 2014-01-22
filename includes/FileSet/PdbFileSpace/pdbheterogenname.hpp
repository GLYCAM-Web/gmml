// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGENNAME_HPP
#define PDBHETEROGENNAME_HPP

#include <string>

namespace PdbFileSpace
{
    class PdbHeterogenName
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbHeterogenName();
            PdbHeterogenName(const std::string& heterogen_identifier, const std::string& heterogen_name);
            PdbHeterogenName(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetHeterogenIdentifier();
            std::string GetHeterogenName();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetHeterogenIdentifier(const std::string heterogen_identifier);
            void SetHeterogenName(const std::string heterogen_name);

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
            std::string heterogen_name_;
    };
}
#endif // PDBHETEROGENNAME_HPP
