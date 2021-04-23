#ifndef SEQUENCE_PARSER_HPP
#define SEQUENCE_PARSER_HPP

#include <string>
#include "includes/InputSet/CondensedSequence/parsedResidue.hpp"

namespace CondensedSequence
{
    class SequenceParser
    {
    public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        SequenceParser(std::string sequence);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::string Print();
        std::vector<ParsedResidue*> GetParsedResidues();
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
    private:
        SequenceParser() {};
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline bool DerivativesExist() {return savedDerivatives_.size();}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        std::vector<std::string> ExtractDerivatives();
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        bool CheckSequenceSanity(std::string sequence);
        void ParseLabelledInput(std::string inString);
        bool ParseCondensedSequence(std::string inString);
        void RecurveParse(size_t &currentIndex, std::string sequence, ParsedResidue* parent);
        ParsedResidue* SaveResidue(const size_t windowStart, const size_t windowEnd, const std::string sequence, ParsedResidue* parent);
        inline void SaveDerivative(std::string derivative) {savedDerivatives_.push_back(derivative);}
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
    	std::vector<std::string> savedDerivatives_;
        std::vector<std::unique_ptr<ParsedResidue>> parsedResidues_;
    };
}
#endif