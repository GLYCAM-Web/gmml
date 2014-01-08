// Author: Alireza Khatamian

#ifndef PDBRESIDUESEQUENCECARD_HPP
#define PDBRESIDUESEQUENCECARD_HPP

#include <string>
#include <map>

namespace PdbFileSpace
{
    class PdbResidueSequence;
    class PdbResidueSequenceCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<char, PdbResidueSequence*> ResidueSequenceMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbResidueSequenceCard();
            PdbResidueSequenceCard(const std::string& record_name);
            PdbResidueSequenceCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            ResidueSequenceMap GetResidueSequenceChain();

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
            ResidueSequenceMap residue_sequence_chains_;

    };
}

#endif // PDBRESIDUESEQUENCECARD_HPP
