#ifndef PDBHETEROGENATOMCARD_HPP
#define PDBHETEROGENATOMCARD_HPP

#include <string>
#include <map>
#include <iostream>

namespace PdbFileSpace
{
    class PdbAtom;

    class PdbHeterogenAtomCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef PdbAtom PdbHeterogenAtom;
            typedef std::map<int, PdbHeterogenAtom*> PdbHeterogenAtomMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbHeterogenAtomCard();
            PdbHeterogenAtomCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            PdbHeterogenAtomMap GetHeterogenAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            PdbHeterogenAtomMap heterogen_atoms_;

    };
}

#endif // PDBHETEROGENATOMCARD_HPP
