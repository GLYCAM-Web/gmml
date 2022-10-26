#include <sstream>
#include <cctype> // isDigit
#include "includes/InputSet/CondensedSequence/assemblyBuilder.hpp"

#include "includes/Abstract/absResidue.hpp" // For the Abstract::ResidueType
#include "includes/MolecularMetadata/GLYCAM/glycam06DerivativeAglyconeInfo.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06DerivativeChargeAdjustment.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06ResidueNameGenerator.hpp" // To get glycam name for ParsedResidue
#include "includes/MolecularModeling/Selections/selections.hpp"
#include "includes/MolecularModeling/assembly.hpp" // Only to use silly Assembly functions. Should go away. 
#include "includes/MolecularModeling/atom.hpp" // For setting Angles and bond distances
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp"

//using Abstract::Residue; // For Abstract::ResidueType

using CondensedSequence::AssemblyBuilder;
using MolecularModeling::Assembly;

AssemblyBuilder::AssemblyBuilder(std::string inputSequence, std::string prepFilePath, Assembly *inputAssembly) : SequenceManipulator{inputSequence} 
{
	this->ReorderSequence(); // Linkages must be in ascending order for looking up Glycam codes? Fix this dependency Oliver.
	this->SetIndexByConnectivity();
    codeUtils::ensureFileExists(prepFilePath);
	PrepFileSpace::PrepFile prepFile(prepFilePath);
	gmml::log(__LINE__,__FILE__,gmml::INF,"Prepfile used is " + prepFilePath);
	this->SetPrepResidueMap(prepFile.GetResidues()); //A mapping between a residue name and its residue object
	this->GenerateResidues(inputAssembly);
	inputAssembly->EnsureIntegralCharge();
	return;
}

void AssemblyBuilder::GenerateResidues(Assembly *assembly)
{
	std::vector<MolecularModeling::Residue> createdResidues;
	createdResidues.reserve(this->getResidues().size());
	CondensedSequence::ParsedResidue *aglycone = this->GetTerminal();
	auto result = this->GetPrepResidueMap()->find(this->GetGlycamResidueName(*aglycone));
	std::stringstream ss;
	ss << "Found prep entry: " << result->first << " for " << aglycone->GetName() << "\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    MolecularModeling::Residue& gmmlParent = assembly->CreateResidue(result->second, aglycone->GetType());
	gmmlParent.addLabel(aglycone->getLabel());
	for (auto &child : aglycone->GetChildren())
	{
		this->RecurveGenerateResidues(child, gmmlParent, assembly);	
	}
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished generating residues for Assembly.\n");
	return;
}

void AssemblyBuilder::RecurveGenerateResidues(ParsedResidue* parsedChild, MolecularModeling::Residue& gmmlParent, 
	Assembly* assembly)
{	
	if (parsedChild->GetType() == Abstract::ResidueType::Deoxy)
	{
		gmml::log(__LINE__, __FILE__, gmml::INF, "Dealing with deoxy for " + gmmlParent.GetName());
		gmmlParent.MakeDeoxy(parsedChild->GetLink());
		return;
	}
	auto prepEntry = this->GetPrepResidueMap()->find(this->GetGlycamResidueName(*parsedChild));
	if (prepEntry == this->GetPrepResidueMap()->end())
	{
		std::string errorMessage = "Could not find prep entry for " + parsedChild->GetName() + ". GLYCAM code used to search is: " + this->GetGlycamResidueName(*parsedChild);
		gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
		throw std::runtime_error(errorMessage);
	}
	else
	{
	    std::stringstream ss;
		ss << "Found prep entry: " << prepEntry->first << " for " << parsedChild->GetName() << "\n";
		gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
		Residue& newGmmlChild = assembly->CreateResidue(prepEntry->second, parsedChild->GetType());
		newGmmlChild.addLabel(parsedChild->getLabel());
		newGmmlChild.SetResidueNumber(parsedChild->getIndex()); // Number and index have become too intertwined. Need to fix this.
		//newGmmlChild.AddEdge(&gmmlParent, parsedChild->GetLinkageName()); // I think this is wishlist versus BondResiduesDeduceAtoms.
		this->BondResiduesDeduceAtoms(gmmlParent, newGmmlChild, parsedChild->GetLinkageName());
		//this->InitializeInterResidueGeometry(gmmlParent, newGmmlChild); // I think this is wishlist versus BondResiduesDeduceAtoms.
		for (auto &child : parsedChild->GetChildren())
		{
			this->RecurveGenerateResidues(child, newGmmlChild, assembly);	
		}
	}
	return;
}

// All this stuff should go into Residue. Residue has private Head and child atoms with private getters/setters. Solves this mess.
// There is way too much checking going on in there. Guarantees need be made elsewhere during construction.
void AssemblyBuilder::BondResiduesDeduceAtoms(MolecularModeling::Residue& parentResidue, MolecularModeling::Residue& childResidue, std::string linkageLabel)
{
	std::stringstream logss;
	logss << "Here with parent " << parentResidue.GetId() << " and child: " << childResidue.GetId() << " and linkageLabel: " << linkageLabel;
	gmml::log(__LINE__,__FILE__,gmml::INF, logss.str());
	// This is using the new Node<Residue> functionality and the old AtomNode 
	parentResidue.addChild(linkageLabel, &childResidue);
	// Now go figure out how which Atoms to bond to each other in the residues.
	// Rule: Can't ever have a child aglycone or a parent derivative.
	std::string parentAtomName, childAtomName;
	if (parentResidue.GetType() == Abstract::ResidueType::Aglycone)
	{ 
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup connectionAtomLookup;
		parentAtomName = connectionAtomLookup.GetConnectionAtomForResidue(parentResidue.GetName());
	}
	else if (parentResidue.GetType() == Abstract::ResidueType::Sugar)
	{ // Linkage example: childb1-4parent, it's never parentb1-4child 
		size_t linkPosition = 3;
		if (childResidue.GetType() == Abstract::ResidueType::Derivative)
		{ // label will be just a single number.
			linkPosition = 0;
		}
		else if (linkageLabel.size() < 4)
		{
		    std::string message = "The deduced linkageLabel is too small:\n" + linkageLabel + ".\nWe require anomer, start atom number, a dash, and connecting atom number. Example:\na1-4";
		    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
		    throw std::runtime_error(message);
		}
		if(!isdigit(linkageLabel.substr(linkPosition).at(0)))
		{
		    std::string message = "Could not convert the last linkage number to an integer: " + linkageLabel;
		    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
		    throw std::runtime_error(message);
		}
		parentAtomName = selection::GetNonCarbonHeavyAtomNumbered(parentResidue.GetAtoms(), linkageLabel.substr(linkPosition));
	}
	else
	{
	    logss << "Error: parent residue: " << parentResidue.GetName() << " with type " << parentResidue.GetType() << " isn't either Aglycone or Sugar, and derivatives cannot be parents.";
	    gmml::log(__LINE__,__FILE__, gmml::ERR, logss.str());
	    throw std::runtime_error(logss.str());
	}
	Atom* parentAtom = parentResidue.GetAtom(parentAtomName);
	if (parentAtom == nullptr)
	{
	    std::string message = "Did not find atom named " + parentAtomName + " in residue: " + parentResidue.GetId();
	    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
	    throw std::runtime_error(message);
	}
	gmml::log(__LINE__,__FILE__,gmml::INF, parentAtom->GetId());
	// Now get child atom
	if (childResidue.GetType() == Abstract::ResidueType::Derivative)
	{
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup connectionAtomLookup;
		childAtomName = connectionAtomLookup.GetConnectionAtomForResidue(childResidue.GetName());	
	}
	else if (childResidue.GetType() == Abstract::ResidueType::Sugar)
	{
		std::string childLinkageNumber = linkageLabel.substr(1,1);
		if(!isdigit(childLinkageNumber.at(0)))
		{
		    std::string message = "Could not convert the first linkage number to an integer: " + childLinkageNumber;
		    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
		    throw std::runtime_error(message);
		}
		childAtomName = "C" + childLinkageNumber;
	}
	else
	{
	    logss << "Error: child residue: " << childResidue.GetName() << " with type " << childResidue.GetType() << " is neither derivative or Sugar (aglycones cannot be children)";
	    gmml::log(__LINE__,__FILE__, gmml::ERR, logss.str());
	    throw std::runtime_error(logss.str());
	}
	Atom* childAtom = childResidue.GetAtom(childAtomName);
	if (childAtom == nullptr)
	{
	    std::string message = "Did not find atom named " + childAtomName + " in residue: " + childResidue.GetId();
	    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
	    throw std::runtime_error(message);
	}
    gmml::log(__LINE__,__FILE__,gmml::INF, childAtom->GetId());
	// Now bond the atoms. Needs to change when AtomNode goes away.
	childAtom->GetNode()->AddNodeNeighbor(parentAtom);
	parentAtom->GetNode()->AddNodeNeighbor(childAtom);
	logss << "Bonded " << parentResidue.GetName() << "@" << parentAtomName << " to " << childResidue.GetName() << "@" << childAtomName << std::endl;
	// Charge adjustment
	if (childResidue.GetType() == Abstract::ResidueType::Derivative)
	{
		logss << "Charge Adjustment.\n";
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeChargeAdjustmentLookupContainer lookup;
		std::string adjustAtomName = lookup.GetAdjustmentAtom(childResidue.GetName());
		adjustAtomName += linkageLabel.substr(0,1); 
		Atom* atomToAdjust = parentResidue.GetAtom(adjustAtomName);
		logss << "    Derivative is " << childResidue.GetName() << ". Adjusting charge on " << atomToAdjust->GetName() << "\n";
		logss << "    Adjusting by: " << lookup.GetAdjustmentCharge(childResidue.GetName()) << "\n";
		gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
		atomToAdjust->SetCharge(atomToAdjust->GetCharge() + lookup.GetAdjustmentCharge(childResidue.GetName()) );
	}
	// Geometry
	logss << "Setting bond distance.\n";
	gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
	MolecularModeling::Assembly whyOhGodWhy; // Doing as few changes as possible. These functions should be in a geometryTopology namespace.	
	whyOhGodWhy.SetResidueResidueBondDistance(parentAtom, childAtom);
	//GeometryTopology::SetDistance(parentAtom, childAtom);
	// Angle
	logss << "Setting angles.\n";
	gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
	const double angle_to_set = 109.4;
	for (auto &parentAtomNeighbor : parentAtom->GetNode()->GetNodeNeighbors())
	{ 
		if ( (parentAtomNeighbor->GetName().at(0) != 'H') && (parentAtomNeighbor != childAtom ) )
		{
			//whyOhGodWhy.SetAngle(parentAtomNeighbor, parentAtom, childAtom, angle_to_set);
			GeometryTopology::SetAngle(parentAtomNeighbor, parentAtom, childAtom, angle_to_set);
		}
	}
	return;	
}
	
std::string AssemblyBuilder::GetGlycamResidueName(ParsedResidue &residue)
{
    std::string linkages = "";
    if (residue.GetType() == Abstract::ResidueType::Sugar)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Checking for glycosidic linkages that connect to " + residue.GetResidueName());
        linkages = residue.GetChildLinkagesForGlycamResidueNaming();
    }
    std::string code = gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNameGenerator(linkages, residue.GetIsomer(), residue.GetResidueName(), 
                                                                            residue.GetRingType(), residue.GetResidueModifier() + residue.GetRingShape(), residue.GetConfiguration() );
    return code;
}
