// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBMODELTYPECARD_HPP
#define PDBMODELTYPECARD_HPP

#include <string>
#include <vector>
#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbModelTypeCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbModelTypeCard();
            PdbModelTypeCard(const std::string& record_name, const std::vector<std::string>& comments);
            PdbModelTypeCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            std::vector<std::string> GetComments();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);
            void SetComments(const std::vector<std::string> comments);

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
            std::string record_name_;
            std::vector<std::string> comments_;
    };
}

#endif // PDBMODELTYPECARD_HPP
