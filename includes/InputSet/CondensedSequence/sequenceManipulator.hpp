#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_SEQUENCE_MANIPULATOR_HPP
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_SEQUENCE_MANIPULATOR_HPP

#include "includes/InputSet/CondensedSequence/parsedResidue.hpp"
#include "includes/InputSet/CondensedSequence/sequenceParser.hpp"
#include "includes/InputSet/CondensedSequence/graphVizDotConfig.hpp"

namespace CondensedSequence
{
    class SequenceManipulator : public SequenceParser
    {
    public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        SequenceManipulator(std::string inputSequence) : SequenceParser{inputSequence} {};
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline ParsedResidue* GetTerminal() {return this->GetParsedResidues().at(0);}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string ReorderSequence();
        void LabelSequence();
        void SetIndexByConnectivity();
        std::string Print(const bool withLabels = false);
        std::vector<ParsedResidue*> GetParsedResiduesOrderedByConnectivity();
        std::string PrintGraphViz(GraphVizDotConfig &configs);
    private:
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void RecurvePrint(ParsedResidue* currentResidue, int& branchStackSize, std::vector<std::string>& output, const bool withLabels);
        std::string GetGraphVizLineForResidue(ParsedResidue &residue, GraphVizDotConfig &configs);
    };
}
#endif
