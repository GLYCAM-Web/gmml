#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_SEQUENCE_MANIPULATOR_HPP
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_SEQUENCE_MANIPULATOR_HPP

#include "includes/InputSet/CondensedSequence/sequenceParser.hpp"

namespace CondensedSequence
{
    struct GraphVizDotConfig 
    {
        public:
            bool show_config_labels_;
            bool show_edge_labels_;
            bool show_position_labels_;
            int dpi_;
            std::string file_name_;
            std::string svg_directory_path_;
            GraphVizDotConfig()
            {
                this->show_edge_labels_ = false;
                this->show_config_labels_ = true;
                this->show_position_labels_ = true;
                this->dpi_ = 72;
                this->svg_directory_path_ = "/programs/gems/gmml/includes/MolecularMetadata/Sugars/SNFG_Symbol_Images/";
                this->file_name_ = "oligosaccharide.dot";
            }
    };
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
        std::string LabelSequence();
        void SetIndexByConnectivity();
        std::string Print(const bool withLabels = false);
        std::string PrintGraphViz(GraphVizDotConfig &configs);
        std::vector<ParsedResidue*> GetParsedResiduesOrderedByConnectivity();
    private:
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void RecurvePrint(ParsedResidue* currentResidue, int& branchStackSize, std::vector<std::string>& output, const bool withLabels);
        std::string GetGraphVizLineForResidue(ParsedResidue &residue, GraphVizDotConfig &configs);
    };
}
#endif
