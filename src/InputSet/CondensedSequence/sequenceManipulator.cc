#include <sstream>
#include "includes/InputSet/CondensedSequence/sequenceManipulator.hpp"
#include "includes/MolecularModeling/Graph/Graph.hpp"

using CondensedSequence::SequenceManipulator;
using CondensedSequence::ParsedResidue;


void SequenceManipulator::ReorderSequence()
{	// Just doing the default by ascending link number for now.
	for (auto &residue : this->GetParsedResidues())
	{
		residue->SortInEdgesBySourceTObjectComparator();
	}
	return;
}

std::vector<ParsedResidue*> SequenceManipulator::GetParsedResiduesOrderedByConnectivity()
{
    std::vector<ParsedResidue*> rawResidues;
    // Go via Graph so order decided by connectivity, depth first traversal:
    TemplateGraph::Graph<ParsedResidue> sequenceGraph(this->GetTerminal());
    for(auto &node : sequenceGraph.GetNodes())
    {
        rawResidues.push_back(node->GetObjectPtr());
    }
    return rawResidues;
}

void SequenceManipulator::LabelSequence()
{
	std::stringstream ss;
	int linkIndex = 0; // Convention to start form 0 for linkages.
	int residueIndex = 1; // Convention to start from 1 for residues.
	//auto startResidue = this->GetTerminal();
	for (auto &residue : this->GetParsedResiduesOrderedByConnectivity())
	{
		ss << residue->GetName() << "&Label=residue-" << residueIndex << ";";
		residue->AddLabel(ss.str());
		++residueIndex;
		ss.str( std::string() ); ss.clear();  // Must do both of these to clear the stream
		for (auto &linkage : residue->GetOutEdges())
		{
			ss << linkage->GetLabel() << "&Label=link-" << linkIndex << ";";
			linkage->AddLabel(ss.str());
			++linkIndex;
			ss.str( std::string() ); ss.clear(); // Must do both of these to clear the stream
		}
	}
	return;
}

// void SequenceManipulator::PrintLabelledSequence()
// {
// 	std::vector<std::string> labelsToPrint;
// 	auto glycamLabelSignature = "&Label=";
// 	for (auto &residue : this->GetParsedResidues())
// 	{
// 		if (residue->GetType() == ParsedResidue::Type::Aglycone)
// 		{ // Aglycone doesn't have linkage, so next for loop doesn't trigger for it.
// 			labelsToPrint.push_back(residue->FindLabelContaining(glycamLabelSignature));
// 		}
// 		for (auto &linkage : residue->GetOutEdges())
// 		{ 
// 			labelsToPrint.push_back(residue->FindLabelContaining(glycamLabelSignature) 
// 				                  + linkage->FindLabelContaining(glycamLabelSignature) );
// 		}
// 	}
// 	std::reverse(labelsToPrint.begin(), labelsToPrint.end()); // Reverse order, as it starts from terminal.
// 	for (auto &label : labelsToPrint)
// 		std::cout << label;
// 	return;
// }

void SequenceManipulator::Print(const bool withLabels)
{
	if (withLabels)
	{
		this->LabelSequence();
	}
	std::vector<std::string> output;
	int branchStackSize = 0;
	this->RecurvePrint(this->GetTerminal(), branchStackSize, output, withLabels);
	std::reverse(output.begin(), output.end()); // Reverse order, as it starts from terminal.
	for (auto &label : output)
		std::cout << label;
}

void SequenceManipulator::RecurvePrint(ParsedResidue* currentResidue, int& branchStackSize, std::vector<std::string>& output, const bool withLabels)
{
	auto neighbors = currentResidue->GetChildren();
	size_t numberOfNeighbors = neighbors.size();
	// Derivatives. E.g. 2S,3Me in DManp[2S,3Me]a1-6DManpa1-OH
	std::string outputResidueString = currentResidue->GetName(withLabels);
	std::vector<std::string> derivatives;
	for (auto &neighbor : neighbors)
	{
		if (neighbor->GetType() == ParsedResidue::Type::Derivative)
		{
			--numberOfNeighbors;
			derivatives.push_back(neighbor->GetLinkageName(withLabels) + neighbor->GetName(withLabels));
			derivatives.push_back(",");
		}
	}
	if (!derivatives.empty())
	{
		derivatives.pop_back(); // Remove the last ","
		outputResidueString += "[";
		for (auto &derivative : derivatives)
		{
			outputResidueString += derivative;
		}
		outputResidueString += "]";
	}
	// Output
	outputResidueString += currentResidue->GetLinkageName(withLabels);
	output.push_back(outputResidueString);
	// End of a branch check
	if (numberOfNeighbors == 0 && branchStackSize > 0)
	{
		output.push_back("["); 
		--branchStackSize;
	}
	size_t loopCount = 0;
	for (auto &neighbor : neighbors)
	{
		if (neighbor->GetType() != ParsedResidue::Type::Derivative)
		{
			++loopCount;
			if (loopCount < numberOfNeighbors)
			{
				output.push_back("]");
				++branchStackSize;
			}
			this->RecurvePrint(neighbor, branchStackSize, output, withLabels);
		}
	}
	return;
}