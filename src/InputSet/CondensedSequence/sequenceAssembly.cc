#include <sstream>
#include "includes/InputSet/CondensedSequence/sequenceAssembly.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06DerivativeAglyconeInfo.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06DerivativeChargeAdjustment.hpp"
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
	std::cout << "\nPREPFILESTUFFFFF\n" << std::endl;
	PrepFileSpace::PrepFile prepFile(prepFilePath);
	//A mapping between a residue name and its residue object
	this->SetPrepResidueMap(prepFile.GetResidues());
	auto aglycone = this->GetTerminal();
	auto result = this->GetPrepResidueMap()->find(aglycone->GetGlycamResidueName());
	std::cout << "Found prep entry: " << result->first << " for " << aglycone->GetName() << "\n";
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

void SequenceAssembly::RecurveGenerateResidues(ParsedResidue* parsedChild, MolecularModeling::Residue& gmmlParent, 
	std::vector<MolecularModeling::Residue>& createdResidues)
{	
	//std::cout << "Recurve Gen Res" << std::endl;
	if (parsedChild->GetType() == ParsedResidue::Type::Deoxy)
	{
		gmmlParent.MakeDeoxy(parsedChild->GetLink());
		return;
	}
	auto prepEntry = this->GetPrepResidueMap()->find(parsedChild->GetGlycamResidueName());
	if (prepEntry == this->GetPrepResidueMap()->end())
	{
		std::cout << "Could not find prep entry for " << parsedChild->GetName() << " Glycam: " << parsedChild->GetGlycamResidueName() << std::endl; 
	}
	else
	{
		std::cout << "Found prep entry: " << prepEntry->first << " for " << parsedChild->GetName() << "\n";
		auto &newGmmlChild = createdResidues.emplace_back(prepEntry->second, parsedChild->GetType());
		newGmmlChild.AddLabel(parsedChild->GetLabel());
		//newGmmlChild.AddEdge(&gmmlParent, parsedChild->GetLinkageName());
		this->BondResiduesDeduceAtoms(gmmlParent, newGmmlChild, parsedChild->GetLinkageName());
		//this->InitializeInterResidueGeometry(gmmlParent, newGmmlChild);
		//std::cout << "Recurve created " << newGmmlChild.GetLabel() << std::endl;
		for (auto &child : parsedChild->GetChildren())
		{
			this->RecurveGenerateResidues(child, newGmmlChild, createdResidues);	
		}
	}
	return;
}

// All this stuff goes to Residue. Residue has private Head and child atoms with private getters/setters. Solves this mess.
void SequenceAssembly::BondResiduesDeduceAtoms(MolecularModeling::Residue& parentResidue, MolecularModeling::Residue& childResidue, std::string linkageLabel)
{
	// This is using the new Node<Residue> functionality and the old AtomNode 
	parentResidue.AddEdge(&childResidue, linkageLabel);
	// Now go figure out how which Atoms to bond to each other in the residues.
	// Rule: Can't ever have a child aglycone or a parent derivative.
	std::string parentAtomName, childAtomName;
	if (parentResidue.GetType() == Residue::Type::Aglycone)
	{ 
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup connectionAtomLookup;
		parentAtomName = connectionAtomLookup.GetConnectionAtomForResidue(parentResidue.GetName());
	}
	else if (parentResidue.GetType() == Residue::Type::Sugar)
	{ // Linkage example: childb1-4parent, it's never parentb1-4child 
		size_t linkPosition = 3;
		if (childResidue.GetType() == Residue::Type::Derivative)
		{ // label will be just a single number.
			linkPosition = 0;
		}
		parentAtomName = selection::GetNonCarbonHeavyAtomNumbered(parentResidue.GetAtoms(), linkageLabel.substr(linkPosition));
	}
	Atom* parentAtom = parentResidue.GetAtom(parentAtomName);
	// Now get child atom
	if (childResidue.GetType() == Residue::Type::Derivative)
	{
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup connectionAtomLookup;
		childAtomName = connectionAtomLookup.GetConnectionAtomForResidue(childResidue.GetName());	
	}
	else if (childResidue.GetType() == Residue::Type::Sugar)
	{
		auto childLinkageNumber = linkageLabel.substr(1,1);
		childAtomName = "C" + childLinkageNumber;
	}
	Atom* childAtom = childResidue.GetAtom(childAtomName);
	// Now bond the atoms. Needs to change when AtomNode goes away.
	childAtom->GetNode()->AddNodeNeighbor(parentAtom);
	parentAtom->GetNode()->AddNodeNeighbor(childAtom);
	std::cout << "Bonded " << parentResidue.GetName() << "@" << parentAtomName << " to " << childResidue.GetName() << "@" << childAtomName << std::endl;
	// Charge adjustment
	if (childResidue.GetType() == Residue::Type::Derivative)
	{
	//	std::cout << "Charge Adjustment.\n";
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeChargeAdjustmentLookupContainer lookup;
		std::string adjustAtomName = lookup.GetAdjustmentAtom(childResidue.GetName());
		adjustAtomName += linkageLabel.substr(0,1); 
		Atom* atomToAdjust = parentResidue.GetAtom(adjustAtomName);
		std::cout << "    Derivative is " << childResidue.GetName() << ". Adjusting charge on " << atomToAdjust->GetName() << "\n";
		std::cout << "    Adjusting by: " << lookup.GetAdjustmentCharge(childResidue.GetName()) << "\n";
		atomToAdjust->SetCharge(atomToAdjust->GetCharge() + lookup.GetAdjustmentCharge(childResidue.GetName()) );
	}
	// Geometry
	//std::cout << "Setting bond distance.\n";
	MolecularModeling::Assembly whyOhGodWhy; // Doing as few changes as possible. These functions should be in a geometryTopology namespace.	
	whyOhGodWhy.SetResidueResidueBondDistance(parentAtom, childAtom);
	// Angle
	//std::cout << "Setting angles.\n";
	const double angle_to_set = 109.4;
	Atom* parentAtomNeighbor;
	for (auto &neighbor : parentAtom->GetNode()->GetNodeNeighbors())
	{ 
		if ( (neighbor->GetName().at(0) != 'H') && (neighbor != childAtom ) )
		{
			parentAtomNeighbor = neighbor;
			//std::cout << "Found neighbor, I'll set the angle now!\n";
			whyOhGodWhy.SetAngle(parentAtomNeighbor, parentAtom, childAtom, angle_to_set);
		}
	}
	return;	
}
	