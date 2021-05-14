#ifndef SEQUENCE_MANIPULATOR_HPP
#define SEQUENCE_MANIPULATOR_HPP

#include "includes/InputSet/CondensedSequence/sequenceParser.hpp"

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
        void ReorderSequence();
        void LabelSequence();
        void SetIndexByConnectivity();
        void PrintLabelledSequence();
        void Print(const bool withLabels = false);
        void PrintGraphViz(std::string SnfgFilePath = "/programs/gems/gmml/includes/MolecularMetadata/Sugars/SNFG_Symbol_Images/");
        std::string GetGraphVizLineForResidue(ParsedResidue &residue, std::string SnfgFilePath);
        std::vector<ParsedResidue*> GetParsedResiduesOrderedByConnectivity();
    private:
        void RecurvePrint(ParsedResidue* currentResidue, int& branchStackSize, std::vector<std::string>& output, const bool withLabels);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
    };
}
#endif