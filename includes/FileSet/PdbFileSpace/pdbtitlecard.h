#ifndef PDBTITLECARD_H
#define PDBTITLECARD_H

#include <string>

namespace PdbFileSpace
{
    class PdbTitleCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbTitleCard();
            PdbTitleCard(std::string record_name, std::string title);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            std::string GetTitle();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(std::string record_name);
            void SetTitle(std::string title);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            std::string title_;
    };
}

#endif // PDBTITLECARD_H
