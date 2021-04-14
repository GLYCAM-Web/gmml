#include <sstream>
#include "includes/InputSet/CondensedSequence/sequenceAssembly.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"

using CondensedSequence::SequenceAssembly;

SequenceAssembly::SequenceAssembly(std::string inputSequence, std::string prepFilePath) : SequenceManipulator{inputSequence} 
{
	this->ReorderSequence(); // Linkages must be in ascending order for looking up Glycam codes? Fix this dependancy Oliver.
	auto molecularModelingResidues = this->GenerateResidues(prepFilePath);
	// std::cout << "Created these molecularModelingResidues: \n";
	// for (auto &mmResidue : molecularModelingResidues)
	// {
	// 	std::cout << mmResidue.GetIndex() << std::endl;
	// }
}

std::vector<MolecularModeling::Residue> SequenceAssembly::GenerateResidues(std::string prepFilePath)
{
	std::vector<MolecularModeling::Residue> createdResidues;
	createdResidues.reserve(this->GetParsedResidues().size());
	std::cout << "\nDoing prepFile stuff" << std::endl;
	PrepFileSpace::PrepFile prepFile(prepFilePath);
	//A mapping between a residue name and its residue object
	this->SetPrepResidueMap(prepFile.GetResidues());
	auto aglycone = this->GetTerminal();
	auto result = this->GetPrepResidueMap()->find(aglycone->GetGlycamResidueName());
	auto &gmmlParent = createdResidues.emplace_back(result->second);
	gmmlParent.AddLabel(aglycone->GetLabel());
	for (auto &child : aglycone->GetChildren())
	{
		this->RecurveGenerateResidues(child, gmmlParent, createdResidues);	
	}
	std::cout << "Created residues: \n";
	for (auto &mmResidue : createdResidues)
	{
		std::cout << mmResidue.GetLabel() << std::endl;
	}
	return createdResidues;
}


void SequenceAssembly::RecurveGenerateResidues( ParsedResidue *currentChild, MolecularModeling::Residue& gmmlParent, 
	std::vector<MolecularModeling::Residue> &createdResidues)
{	
	auto prepEntry = this->GetPrepResidueMap()->find(currentChild->GetGlycamResidueName());
	if (prepEntry == this->GetPrepResidueMap()->end())
	{
		std::cout << "Could not find prep entry.\n"; 
	}
	else
	{
		auto &newGmmlChild = createdResidues.emplace_back(prepEntry->second);
		newGmmlChild.AddLabel(currentChild->GetLabel());
		newGmmlChild.AddEdge(&gmmlParent, currentChild->GetLinkageLabel());
		std::cout << "Recurve created " << newGmmlChild.GetLabel() << std::endl;
		for (auto &child : currentChild->GetChildren())
		{
			this->RecurveGenerateResidues(child, newGmmlChild, createdResidues);	
		}
	}
	return;
}

// MolecularModeling::Residue* SequenceAssembly::CreateResidue(PrepFileSpace::PrepFileResidue *prep_residue, MolecularModeling::Residue parent, std::string linkageLabel)
// {
// 	auto residues = this->GetNewResidues();
// 	auto newChildResidue = residues.emplace_back(prep_residue);
// 	newChildResidue.AddLinkage(parent, linkageLabel);
// 	return newChildResidue;
// }