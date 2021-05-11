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
	this->SetIndexByConnectivity();
	std::stringstream ss;
	for (auto &residue : this->GetParsedResiduesOrderedByConnectivity())
	{
		ss << residue->GetName() << "&Label=residue-" << residue->GetIndex() << ";";
		residue->AddLabel(ss.str());
		ss.str( std::string() ); ss.clear();  // Must do both of these to clear the stream
		for (auto &linkage : residue->GetOutEdges())
		{
			ss << linkage->GetLabel() << "&Label=link-" << linkage->GetIndex() << ";";
			linkage->AddLabel(ss.str());
			ss.str( std::string() ); ss.clear(); // Must do both of these to clear the stream
		}
	}
	return;
}

void SequenceManipulator::SetIndexByConnectivity()
{
	unsigned long long linkIndex = 0; // Convention to start form 0 for linkages.
	unsigned long long residueIndex = 1; // Convention to start from 1 for residues.
	for (auto &residue : this->GetParsedResiduesOrderedByConnectivity())
	{
		residue->SetIndex(residueIndex);
		++residueIndex;
		for (auto &linkage : residue->GetOutEdges())
		{
			linkage->SetIndex(linkIndex);
			++linkIndex;
		}
	}
	return;
}

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

void SequenceManipulator::PrintGraphViz(std::string SnfgFilePath)
{
	this->SetIndexByConnectivity();
	std::string output = "graph G {graph [splines=false forcelabels=true  dpi=72];\n";
	output += "node [ shape=\"none\" fontname=Helvetica labelfontsize=12 forcelabels=\"true\";\n";
	output += "label=\"none\" size=50 fixedsize=\"true\" scale=\"true\"];\n";
	output += "edge [labelfontsize=12 fontname=Helvetica labeldistance=1.2 labelangle = 320.0];\n";
	output += "rankdir=LR nodesep=\"0.05\" ranksep=\"0.8\";\n";
	for (auto &residue : this->GetParsedResiduesOrderedByConnectivity())
	{
		if (residue->GetType() != ParsedResidue::Type::Derivative) 
		{
			output += residue->GetGraphVizLine(SnfgFilePath) + "\n";
		}
	}
	output += "}\n";
	std::cout << output;
}