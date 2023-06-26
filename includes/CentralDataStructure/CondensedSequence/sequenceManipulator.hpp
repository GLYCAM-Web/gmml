#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphVizDotConfig.hpp"

namespace cdsCondensedSequence
{
    class SequenceManipulator : public SequenceParser
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        SequenceManipulator(std::string inputSequence) : SequenceParser {inputSequence} {};

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline ParsedResidue* GetTerminal() const
        {
            return static_cast<ParsedResidue*>(this->getResidues().front());
        }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string ReorderSequence();
        void LabelSequence();
        void SetIndexByConnectivity();
        std::string Print(const bool withLabels = false) const;
        std::vector<ParsedResidue*> GetParsedResiduesOrderedByConnectivity() const;
        std::string PrintGraphViz(GraphVizDotConfig& configs);

      private:
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void RecurvePrint(ParsedResidue* currentResidue, int& branchStackSize, std::vector<std::string>& output,
                          const bool withLabels) const;
        std::string GetGraphVizLineForResidue(ParsedResidue& residue, GraphVizDotConfig& configs);
    };
} // namespace cdsCondensedSequence
#endif
