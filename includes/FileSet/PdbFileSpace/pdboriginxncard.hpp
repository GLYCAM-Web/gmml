#ifndef PDBORIGINXNCARD_HPP
#define PDBORIGINXNCARD_HPP

#include <string>
#include <vector>
#include <sstream>

namespace PdbFileSpace
{
    class PdbOriginXn;

    class PdbOriginXnCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::vector< PdbOriginXn* > OriginXnVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbOriginXnCard();
            PdbOriginXnCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            OriginXnVector GetOriginXN();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetOriginXN(const OriginXnVector origin_x_n);
            void AddOriginXN(PdbOriginXn* origin);

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
            OriginXnVector origin_x_n_;

    };
}

#endif // PDBORIGINXNCARD_HPP
