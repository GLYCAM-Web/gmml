#include "includes/InputSet/CondensedSequence/carbohydrate.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06DerivativeAglyconeInfo.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06DerivativeChargeAdjustment.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06ResidueNameGenerator.hpp" // To get glycam name for ParsedResidue
#include "includes/Abstract/absResidue.hpp" // For the Residue::Type
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/ParameterSet/PrepFile/prepFile.hpp"
#include "includes/ParameterSet/PrepFile/prepResidue.hpp"
#include "includes/ParameterSet/PrepFile/prepAtom.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Shapers/geometryTopologyInterface.hpp"
#include <sstream>
#include <cctype> // isDigit
//using Abstract::absResidue; // For Residue::Type
using CondensedSequence::Carbohydrate;

Carbohydrate::Carbohydrate(std::string inputSequence, std::string prepFilePath) : SequenceManipulator{inputSequence}
{
	this->ReorderSequence(); // Linkages must be in ascending order for looking up Glycam codes? Fix this dependency Oliver. Update: Fixed. Todo: Confirm/Test the fix Oliver. Just delete this line and test you toolbag.
	this->SetIndexByConnectivity(); // For reporting residue index numbers to the user
	prep::PrepFile glycamPrepFileSelect(prepFilePath, this->GetGlycamNamesOfResidues()); // PrepFile is a cds::Molecule i.e. contains a vector of residues.
	for( auto &cdsResidue: this->getResidues() )
	{
	    ParsedResidue* parsedResidue = static_cast<ParsedResidue*>(cdsResidue);
	    std::cout << "parsedResidue is " << parsedResidue->getName() << std::endl;
	    std::cout << "it's glycam name is " << this->GetGlycamResidueName(parsedResidue) << std::endl;
        std::cout << "it's label is " << parsedResidue->getLabel() << std::endl;
        std::cout << parsedResidue->PrintToString() << std::endl;
        // Deoxy
        if (parsedResidue->GetType() == Abstract::ResidueType::Deoxy)
        {
            std::cout << "Dealing with deoxy for " << parsedResidue->getName() << "\n";
            gmml::log(__LINE__, __FILE__, gmml::INF, "Dealing with deoxy for " + parsedResidue->getName());
            cds::Residue* residueToBeDeoxified = parsedResidue->getParents().at(0);
            residueToBeDeoxified->MakeDeoxy(parsedResidue->GetLink());
            this->deleteResidue(parsedResidue); // Remove the deoxy derivative now.
        }
        else // Not Deoxy
        {
            cds::Residue* prepResidue = glycamPrepFileSelect.getResidue(this->GetGlycamResidueName(parsedResidue));
            if (prepResidue == nullptr)
            {
                std::string message = "Did not find prep entry for " + parsedResidue->getName() + " with glycam residue code: " + this->GetGlycamResidueName(parsedResidue);
                gmml::log(__LINE__,__FILE__,gmml::ERR, message);
                throw(std::runtime_error(message));
            }
            std::cout << "prepResidue is " << prepResidue->getName() << std::endl;
            // ParsedResidue has the correct residue connectivities, prepResidue has the atoms.
            parsedResidue->setAtoms(prepResidue->extractAtoms()); // This moves the atoms, i.e. for prepResidue "Moved from objects are left in a valid but unspecified state"
            glycamPrepFileSelect.deleteResidue(prepResidue); // Death to the prepResidue, if there are repeats with the same name, the next search would find the one without atoms.
            std::cout << "Finished moving atoms from prepResidue to parsed Residue. Adventure awaits! Huzzah!" << std::endl;
            // Charge adjustment. If residue is Derivaitve, the parent residue connection atom's charge is adjusted.
            if (parsedResidue->GetType() == Abstract::ResidueType::Derivative)
            {
                std::stringstream ss;
                ss << "Charge Adjustment.\n";
                gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeChargeAdjustmentLookupContainer lookup;
                std::string adjustAtomName = lookup.GetAdjustmentAtom(parsedResidue->getName());
                adjustAtomName += parsedResidue->GetLinkageName().substr(0,1);
                cds::Atom* atomToAdjust = parsedResidue->getParents().at(0)->FindAtom(adjustAtomName);
                ss << "    Derivative is " << parsedResidue->getName() << ". Adjusting charge on " << atomToAdjust->getName() << "\n";
                ss << "    Adjusting by: " << lookup.GetAdjustmentCharge(parsedResidue->getName()) << "\n";
                gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
                atomToAdjust->setCharge(atomToAdjust->getCharge() + lookup.GetAdjustmentCharge(parsedResidue->getName()) );
            }
        }
	}
    std::cout << "\n\n\nOn to setting 3d structure!\n\n";
	for( auto &cdsResidue: this->getResidues() )
	{
	    for( auto &parentNeighbor : cdsResidue->getParents()) //
	    {
	        std::cout << "Setting connection between " << cdsResidue->getName() << " and it's parent " << parentNeighbor->getName() << ", which has linkageLabel: "
	                << static_cast<ParsedResidue*>(parentNeighbor)->GetLinkageName() << "\n";
	        this->BondResiduesDeduceAtoms(cdsResidue, parentNeighbor);
	    }
	}
	std::cout << "\n\nAre there any premature dtors here today?\n\n";
	//Ensure integralCharge can be a free function that accepts atom vector right?
//	this->EnsureIntegralCharge(inputAssembly->GetTotalCharge());
	return;
}

void Carbohydrate::EnsureIntegralCharge(double charge)
{
	std::stringstream ss;
	ss << std::fixed;
	ss << "Total charge is: " << std::setprecision(5) << charge << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
	double difference = std::fabs(charge - (std::round(charge)));
	if (difference > 0.00001 && difference < 0.99999)
	{
		std::stringstream errorMessage; 
		errorMessage << "Non-integral charge (" << charge << "). You cannot run MD with this.\n";
		std::cerr << errorMessage.str();
		throw errorMessage.str();
	}
	return;
}

//void Carbohydrate::GenerateResidues(Assembly *assembly)
//{
//	std::vector<MolecularModeling::Residue> createdResidues;
//	createdResidues.reserve(this->GetParsedResidues().size());
//	CondensedSequence::ParsedResidue *aglycone = this->GetTerminal();
//	auto result = this->GetPrepResidueMap()->find(this->GetGlycamResidueName(*aglycone));
//	std::stringstream ss;
//	ss << "Found prep entry: " << result->first << " for " << aglycone->getName() << "\n";
//    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
//	Residue& gmmlParent = assembly->CreateResidue(result->second, aglycone->GetType());
//	gmmlParent.addLabel(aglycone->getLabel());
//	for (auto &child : aglycone->GetChildren())
//	{
//		this->RecurveGenerateResidues(child, gmmlParent, assembly);
//	}
//    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished generating residues for Assembly.\n");
//	return;
//}
//
//void Carbohydrate::RecurveGenerateResidues(ParsedResidue* parsedChild, MolecularModeling::Residue& gmmlParent,
//	Assembly* assembly)
//{
//	if (parsedChild->GetType() == ParsedResidue::Type::Deoxy)
//	{
//		gmml::log(__LINE__, __FILE__, gmml::INF, "Dealing with deoxy for " + gmmlParent.getName());
//		gmmlParent.MakeDeoxy(parsedChild->GetLink());
//		return;
//	}
//	auto prepEntry = this->GetPrepResidueMap()->find(this->GetGlycamResidueName(*parsedChild));
//	if (prepEntry == this->GetPrepResidueMap()->end())
//	{
//		std::string errorMessage = "Could not find prep entry for " + parsedChild->getName() + ". GLYCAM code used to search is: " + this->GetGlycamResidueName(*parsedChild);
//		gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
//		throw std::runtime_error(errorMessage);
//	}
//	else
//	{
//	    std::stringstream ss;
//		ss << "Found prep entry: " << prepEntry->first << " for " << parsedChild->getName() << "\n";
//		gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
//		Residue& newGmmlChild = assembly->CreateResidue(prepEntry->second, parsedChild->GetType());
//		newGmmlChild.addLabel(parsedChild->getLabel());
//		newGmmlChild.SetResidueNumber(parsedChild->getIndex()); // Number and index have become too intertwined. Need to fix this.
//		//newGmmlChild.AddEdge(&gmmlParent, parsedChild->GetLinkageName()); // I think this is wishlist versus BondResiduesDeduceAtoms.
//		this->BondResiduesDeduceAtoms(gmmlParent, newGmmlChild, parsedChild->GetLinkageName());
//		//this->InitializeInterResidueGeometry(gmmlParent, newGmmlChild); // I think this is wishlist versus BondResiduesDeduceAtoms.
//		for (auto &child : parsedChild->GetChildren())
//		{
//			this->RecurveGenerateResidues(child, newGmmlChild, assembly);
//		}
//	}
//	return;
//}

void Carbohydrate::BondResiduesDeduceAtoms(cds::Residue* childResidue, cds::Residue* parentResidue)
{
    using Abstract::ResidueType;
    using cds::Atom;
    std::string linkageLabel = static_cast<ParsedResidue*>(childResidue)->GetLinkageName();
	gmml::log(__LINE__,__FILE__,gmml::INF, "Here with child " + childResidue->getId() + " and parent: " + parentResidue->getId() + " and linkageLabel: " + linkageLabel);
	// This is using the new Node<Residue> functionality and the old AtomNode
	// Now go figure out how which Atoms to bond to each other in the residues.
	// Rule: Can't ever have a child aglycone or a parent derivative.
	Atom* parentAtom = nullptr;
	std::string childAtomName, parentAtomName;
	if (parentResidue->GetType() == ResidueType::Aglycone)
	{
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup connectionAtomLookup;
		parentAtomName = connectionAtomLookup.GetConnectionAtomForResidue(parentResidue->getName());
		parentAtom = parentResidue->FindAtom(parentAtomName);
	}
	else if (parentResidue->GetType() == ResidueType::Sugar)
	{ // Linkage example: childb1-4parent, it's never parentb1-4child
		size_t linkPosition = 3;
		if (childResidue->GetType() == ResidueType::Derivative)
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
		parentAtom = TemplatedSelections::getNonCarbonHeavyAtomNumbered(parentResidue->getAtoms(), linkageLabel.substr(linkPosition));
		std::cout << "Parent atom is " << parentAtom->getName() << std::endl;
	}
	else
	{
	    std::string message = "Error: parent residue: " + parentResidue->getName() + " isn't either Aglycone or Sugar, and derivatives cannot be parents.";
	    gmml::log(__LINE__,__FILE__, gmml::ERR, message);
	    throw std::runtime_error(message);
	}
	if (parentAtom == nullptr)
	{
	    std::string message = "Did not find connection atom in residue: " + parentResidue->getId();
	    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
	    throw std::runtime_error(message);
	}
	gmml::log(__LINE__,__FILE__,gmml::INF, parentAtom->getId());
	// Now get child atom
	if (childResidue->GetType() == ResidueType::Derivative)
	{
		gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup connectionAtomLookup;
		std::string glycamNameForResidue = this->GetGlycamResidueName(static_cast<ParsedResidue*>(childResidue));
		childAtomName = connectionAtomLookup.GetConnectionAtomForResidue(glycamNameForResidue);
	}
	else if (childResidue->GetType() == ResidueType::Sugar)
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
	    std::string message = "Error: child residue: " + childResidue->getName() + " is neither derivative or Sugar (aglycones cannot be children)";
	    gmml::log(__LINE__,__FILE__, gmml::ERR, message);
	    throw std::runtime_error(message);
	}
	Atom* childAtom = childResidue->FindAtom(childAtomName);
	if (childAtom == nullptr)
	{
	    std::string message = "Did not find child atom named " + childAtomName + " in child residue: " + childResidue->getId();
	    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
	    throw std::runtime_error(message);
	}
    gmml::log(__LINE__,__FILE__,gmml::INF, childAtom->getId());
    std::stringstream logss;
	// Geometry
	logss << "Setting bond distance between parent " << parentAtom->getName() << " and child " << childAtom->getName() << ".\n";
	gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
	GeometryTopology::FindAtomsToMoveSetDistance(parentAtom, childAtom);
	//   Now bond the atoms. This could also set distance?, and angle? if passed to function?
	childAtom->addBond(parentAtom); // parentAtom also connected to childAtom. Fancy.
	logss << "Bonded " << parentResidue->getName() << "@" << parentAtom->getName() << " to " << childResidue->getName() << "@" << childAtomName << std::endl;
	// Angle
	logss << "Setting angles.\n";
	gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
	const double angle_to_set = 109.4;
	for (auto &parentAtomNeighbor : parentAtom->getNeighbors())
	{
		if ( (parentAtomNeighbor->getName().at(0) != 'H') && (parentAtomNeighbor != childAtom ) )
		{
		   // std::cout << "Setting angle between parentNeighbor " << parentAtomNeighbor->getName() << ", parent " << parentAtom->getName() << " and child " << childAtom->getName() << ". child residue " << childResidue->getName() << " will move\n";
			GeometryTopology::SetAngle(parentAtomNeighbor->getCoordinate(), parentAtom->getCoordinate(), childAtom->getCoordinate(), angle_to_set, childResidue->getCoordinates());
			break;
		}
	}

	return;
}

std::vector<std::string> Carbohydrate::GetGlycamNamesOfResidues() const
{
	std::vector<std::string> names(this->getResidues().size()); // set size of vec for speed.
    std::cout << "Glycam names are: ";
	for(auto &residue : this->getResidues())
	{
	    if (residue->GetType() != Abstract::ResidueType::Deoxy)
	    {
	        names.push_back(this->GetGlycamResidueName(static_cast<ParsedResidue*>(residue)));
	        std::cout << names.back() << ", ";
	    }
	}
	std::cout << "\n";
	return names;
}


std::string Carbohydrate::GetGlycamResidueName(ParsedResidue *residue) const
{
    if (residue->GetType() == Abstract::ResidueType::Deoxy)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Bad idea: We asked for Glycam Residue Name of a deoxy type residuw with name: " + residue->GetResidueName());
        return "";
    }
    std::string linkages = "";
    if (residue->GetType() == Abstract::ResidueType::Sugar)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Checking for glycosidic linkages that connect to " + residue->GetResidueName());
        linkages = residue->GetChildLinkagesForGlycamResidueNaming();
    }
    std::string code = gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNameGenerator(linkages, residue->GetIsomer(), residue->GetResidueName(),
                                                                            residue->GetRingType(), residue->GetResidueModifier() + residue->GetRingShape(), residue->GetConfiguration() );
    return code;
}
