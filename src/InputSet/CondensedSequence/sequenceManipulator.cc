#include <sstream>
#include "includes/InputSet/CondensedSequence/sequenceManipulator.hpp"

using CondensedSequence::SequenceManipulator;
using CondensedSequence::ParsedResidue;


void SequenceManipulator::ReorderSequence()
{
	// Just doing the default by ascending link number for now.
	for (auto &residue : this->GetParsedResidues())
	{
		std::cout << "Links for " << residue->GetName() << ":\n";
		residue->SortInEdgesBySourceTObjectComparator();
		for (auto &neighbor : residue->GetChildren())
		{
			std::cout << neighbor->GetLink() << ",";
		}
		std::cout << "\n";
	}
	return;
}

void SequenceManipulator::LabelSequence()
{
	std::stringstream ss;
	int linkIndex = 0; // Convention to start form 0 for linkages.
	int residueIndex = 1; // Convention to start from 1 for residues.
	auto startResidue = this->GetTerminal();
	for (auto &residue : this->GetParsedResiduesOrderedByConnectivity())
	{
		ss << residue->GetName() << "&Label=residue-" << residueIndex << ";";
		residue->AddLabel(ss.str());
		++residueIndex;
		ss.str( std::string() ); ss.clear();  // Must do both of these to clear the stream
		for (auto &linkage : residue->GetOutEdges())
		{
			std::cout << linkage->GetLabel() << std::endl;
			ss << linkage->GetLabel() << "&Label=link-" << linkIndex << ";";
			linkage->AddLabel(ss.str());
			++linkIndex;
			ss.str( std::string() ); ss.clear(); // Must do both of these to clear the stream
		}
	}
	return;
}

void SequenceManipulator::PrintLabelledSequence()
{
	std::vector<std::string> labelsToPrint;
	auto glycamLabelSignature = "&Label=";
	for (auto &residue : this->GetParsedResidues())
	{
		if (residue->GetType() == ParsedResidue::Type::Aglycone)
		{ // Aglycone doesn't have linkage, so next for loop doesn't trigger for it.
			labelsToPrint.push_back(residue->FindLabelContaining(glycamLabelSignature));
		}
		for (auto &linkage : residue->GetOutEdges())
		{ 
			labelsToPrint.push_back(residue->FindLabelContaining(glycamLabelSignature) 
				                  + linkage->FindLabelContaining(glycamLabelSignature) );
		}
	}
	std::reverse(labelsToPrint.begin(), labelsToPrint.end()); // Reverse order, as it starts from terminal.
	for (auto &label : labelsToPrint)
		std::cout << label;
	return;
}