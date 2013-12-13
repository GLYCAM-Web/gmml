#ifndef PDBHEADERCARD_H
#define PDBHEADERCARD_H

#include <string>

namespace PdbFileSpace
{
    class PdbHeaderCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbHeaderCard();
            PdbHeaderCard(const std::string& record_name, const std::string& classification, const std::string& deposition_date, const std::string& identifier_code);

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
            void SetRecordName(const std::string record_name);
            void SetClassification(const std::string classification);
            void SetDepositionDate(const std::string deposition_date);
            void SetIdentificationCode(const std::string identifier_code);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

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
