#ifndef PDBHEADERCARD_H
#define PDBHEADERCARD_H

#include <string>

namespace PdbFileSpace
{
    class PdbHeaderClass
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PdbHeaderClass();
            PdbHeaderClass(std::string record_name, std::string classification, std::string deposition_date, std::string identifier_code);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            std::string GetClassification();
            std::string GetDepositionDate();
            std::string GetIdentifierCode();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(std::string record_name);
            void SetClassification(std::string classification);
            void SetDepositionDate(std::string deposition_date);
            void SetIdentificationCode(std::string identifier_code);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////



        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            std::string classification_;
            std::string deposition_date_;
            std::string identifier_code_;

    };
}
#endif // PDBHEADERCARD_H
