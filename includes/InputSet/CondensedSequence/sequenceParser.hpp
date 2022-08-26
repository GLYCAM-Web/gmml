#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_SEQUENCE_PARSER_HPP
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_SEQUENCE_PARSER_HPP

#include <string>
#include "includes/InputSet/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/cdsMolecule.hpp"
#include "includes/CentralDataStructure/atom.hpp"

namespace CondensedSequence
{
    class SequenceParser : public cds::cdsMolecule<ParsedResidue, cds::Atom>
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
        std::string parseRepeatingUnits(const std::string inputSequence);
        size_t seekRepeatStart(const std::string &inputSequence, size_t startPosition);
        bool ParseCondensedSequence(std::string inString);
        void RecurveParse(size_t &currentIndex, std::string sequence, ParsedResidue* parent);
        void RecurveParseAlt(size_t &currentIndex, std::string sequence, ParsedResidue* parent);
        ParsedResidue* SaveResidue(const size_t windowStart, const size_t windowEnd, const std::string sequence, ParsedResidue* parent);
        inline void SaveDerivative(std::string derivative) {savedDerivatives_.push_back(derivative);}
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
    	std::vector<std::string> savedDerivatives_;
    };
}
#endif
