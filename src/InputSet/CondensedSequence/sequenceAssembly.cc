#include <sstream>
#include "includes/InputSet/CondensedSequence/sequenceAssembly.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06DerivativeAglyconeInfo.hpp"
#include "includes/MolecularModeling/Abstract/Residue.hpp" // For the Residue::Type
#include "includes/MolecularModeling/Selections/selections.hpp"
#include "includes/MolecularModeling/assembly.hpp" // Only to use silly Assembly functions. Should go away. 
#include "includes/MolecularModeling/atom.hpp" // For setting Angles and bond distances

//using Abstract::Residue; // For Residue::Type

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
	auto &gmmlParent = createdResidues.emplace_back(result->second, aglycone->GetType());
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

void SequenceAssembly::RecurveGenerateResidues( ParsedResidue *parsedChild, MolecularModeling::Residue& gmmlParent, 
	std::vector<MolecularModeling::Residue> &createdResidues)
{	
	auto prepEntry = this->GetPrepResidueMap()->find(parsedChild->GetGlycamResidueName());
	if (prepEntry == this->GetPrepResidueMap()->end())
	{
		std::cout << "Could not find prep entry.\n"; 
	}
	else
	{
		auto &newGmmlChild = createdResidues.emplace_back(prepEntry->second, parsedChild->GetType());
		newGmmlChild.AddLabel(parsedChild->GetLabel());
		//newGmmlChild.AddEdge(&gmmlParent, parsedChild->GetLinkageLabel());
		this->BondResiduesDeduceAtoms(gmmlParent, newGmmlChild, parsedChild->GetLinkageLabel());
		this->InitializeInterResidueGeometry(gmmlParent, newGmmlChild);
		std::cout << "Recurve created " << newGmmlChild.GetLabel() << std::endl;
		for (auto &child : parsedChild->GetChildren())
		{
			this->RecurveGenerateResidues(child, newGmmlChild, createdResidues);	
		}
	}
	return;
}

void SequenceAssembly::BondResiduesDeduceAtoms(MolecularModeling::Residue& parent, MolecularModeling::Residue& child, std::string linkageLabel)
{
	// This is using the new Node<Residue> functionality and the old AtomNode 
	parent.AddEdge(&child, linkageLabel);
	std::cout << "Have entered with types " << parent.GetType() << ", " << child.GetType() << "\n";
	std::cout << "Linkage label: " << linkageLabel << "\n";
	// Now go figure out how which Atoms to bond to each other in the residues.
	// Rule: Can't ever have a child aglycone or a parent derivative.
	std::string parentAtomName, childAtomName;
	if (parent.GetType() == Residue::Type::Aglycone)
	{ 
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup connectionAtomLookup;
		parentAtomName = connectionAtomLookup.GetConnectionAtomForResidue(parent.GetName());
	}
	else if (parent.GetType() == Residue::Type::Sugar)
	{ // Linkage example: child1-4parent. It's never parent1-4child.
		size_t linkPosition = 2;
		if (child.GetType() == Residue::Type::Derivative)
		{ // label will be just a single number.
			linkPosition = 0;
		}
		auto parentLinkageNumber = linkageLabel.at(linkPosition);
		parentAtomName = selection::GetNonCarbonHeavyAtomNumbered(parent.GetAtoms(), parentLinkageNumber);
	}
	std::cout << "Found parent atom " << parentAtomName << "\n";
	auto parentAtom = parent.GetAtom(parentAtomName);

	// Now get child atom
	if (child.GetType() == Residue::Type::Derivative)
	{
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup connectionAtomLookup;
		childAtomName = connectionAtomLookup.GetConnectionAtomForResidue(child.GetName());	
	}
	else if (child.GetType() == Residue::Type::Sugar)
	{
		auto childLinkageNumber = linkageLabel.substr(0,1);
		childAtomName = "C" + childLinkageNumber;
	}
	std::cout << "Found child atom " << childAtomName << "\n";
	auto childAtom = child.GetAtom(childAtomName);
	// Now bond the atoms. Needs to change when AtomNode goes away.
	childAtom->GetNode()->AddNodeNeighbor(parentAtom);
	parentAtom->GetNode()->AddNodeNeighbor(childAtom);
	std::cout << "Bonded " << parent.GetName() << "@" << parentAtomName << " to "
						   << child.GetName() << "@" << childAtomName << std::endl;
	// Now that I have these atoms, I'm going to do geometry stuff
	MolecularModeling::Assembly whyOhGodWhy; // Doing as few changes as possible. These functions should be in a geometryTopology namespace.	
	whyOhGodWhy.SetResidueResidueBondDistance(parentAtom, childAtom);
	// Angle
	Atom* parentAtomNeighbor;
	for (auto &neighbor : parentAtom->GetNode()->GetNodeNeighbors())
	{ 
		if ( (neighbor->GetName().at(0) != 'H') && (neighbor != childAtom ) )
		{
			parentAtomNeighbor = neighbor;
		}
	}
	const double angle_to_set = 109.4;
    whyOhGodWhy.SetAngle(parentAtomNeighbor, parentAtom, childAtom, angle_to_set);
	return;	
}
	
	// Bond that to anomericCarbon in child. Use AtomNode for now. Change later when Node<Atom>.

void SequenceAssembly::InitializeInterResidueGeometry(MolecularModeling::Residue& parent, MolecularModeling::Residue& child)
{

	// 
}