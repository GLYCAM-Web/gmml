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
#include "includes/CentralDataStructure/Writers/cdsOffWriter.hpp"
#include <sstream>
#include <cctype> // isDigit
//using Abstract::absResidue; // For Residue::Type
using CondensedSequence::Carbohydrate;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
// Note: SequenceManipulator constructor can throw. Hmm.
Carbohydrate::Carbohydrate(std::string inputSequence, std::string prepFilePath) : SequenceManipulator{inputSequence}
{
    try
    {
        // Better to throw once I figure out how to catch it in gems. This setting status thing and checking it is a bad pattern.
        this->ReorderSequence(); // Linkages must be in ascending order for looking up Glycam codes? Fix this dependency Oliver. Update: Fixed. Todo: Confirm/Test the fix Oliver. Just delete this line and test you toolbag.
        this->SetIndexByConnectivity(); // For reporting residue index numbers to the user
        // Find relevant Prep residues:
        prep::PrepFile glycamPrepFileSelect(prepFilePath, this->GetGlycamNamesOfResidues()); // PrepFile is a cds::Molecule i.e. contains a vector of residues.
        for( auto &cdsResidue: this->getResidues() )
        {
            // Move atoms from prep file into parsedResidues.
            if (cdsResidue->GetType() != Abstract::ResidueType::Deoxy)
            {
                ParsedResidue* parsedResidue = static_cast<ParsedResidue*>(cdsResidue);
                this->MoveAtomsFromPrepResidueToParsedResidue(glycamPrepFileSelect, parsedResidue);
                // Deal with adjusting charges for derivatives
                if (parsedResidue->GetType() == Abstract::ResidueType::Derivative)
                {
                    this->DerivativeChargeAdjustment(parsedResidue);
                }
            }
        }
        // Apply any deoxy
        for( auto &cdsResidue: this->getResidues() )
        {
            if( cdsResidue->GetType() == Abstract::ResidueType::Deoxy)
            {
                this->ApplyDeoxy(static_cast<ParsedResidue*>(cdsResidue));
            }
        }
        std::cout << "\n\n\nOn to setting 3d structure!\n\n";
        // Set 3D structure
        for( auto &cdsResidue: this->getResidues() )
        {
            for( auto &parentNeighbor : cdsResidue->getParents()) //
            {
                std::cout << "Setting connection between " << cdsResidue->getName() << " and it's parent " << parentNeighbor->getName() << ", which has linkageLabel: "
                        << static_cast<ParsedResidue*>(parentNeighbor)->GetLinkageName() << "\n";
                this->ConnectAndSetGeometry(cdsResidue, parentNeighbor);
            }
        }
        // Set torsions now that everything is attached and the molecule is whole.
//        for( auto &cdsResidue: this->getResidues() )
//        {
//            for( auto &parentNeighbor : cdsResidue->getParents()) //
//            {
//                std::cout << "Finding rotatable dihedrals and applying metadata." << std::endl;
//                cds::ResidueLinkage& linkage = glycosidicLinkages_.emplace_back(cdsResidue, parentNeighbor);
//                std::cout << "Setting default shape" << std::endl;
//                linkage.SetDefaultShapeUsingMetadata();
//            }
//        }
//        this->Generate3DStructureFiles("./", "defaultGeometry");
//        // Wiggle to resolve overlaps:
//        std::cout << "Resolving overlaps" << std::endl;
//        std::vector<cds::Atom*> allAtomsInCarb = this->getAtoms();
//        this->ResolveOverlaps();
//        std::cout << "Overlaps resolved" << std::endl;

        // Ok if have done greedy then the atoms to move needs to beupdated for every linkage:
        std::cout << "Re-determining atoms that need to move for each linkage:" << std::endl;
        for (auto &linkage : glycosidicLinkages_)
        {
            linkage.DetermineAtomsThatMove();
        }
        std::cout << "Final overlap resolution" << std::endl;
        this->ResolveOverlaps();
        std::cout << "Overlaps resolved" << std::endl;
        std::cout << "Number of residues is " << this->getResidues().size() << "\n";
    }
    catch(const std::string &exceptionMessage)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "carbohydrateBuilder class caught this exception message: " + exceptionMessage);
        this->SetStatus("ERROR", exceptionMessage);
        return;
    }
    catch (const std::runtime_error &error)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
        this->SetStatus("ERROR", error.what());
        return;
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "carbohydrateBuilder class caught a throw that was not anticipated. Curious. Death cometh?");
        this->SetStatus("ERROR", "carbohydrateBuilder caught a throw type that was not anticipated. Pretty please report how you got to this to glycam@gmail.com.");
        return;
    }
    //	Ensure integralCharge can be a free function that accepts atom vector right?
    //	this->EnsureIntegralCharge(inputAssembly->GetTotalCharge());
    return;
}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
//std::string fileOutputDirectory = "unspecified", std::string fileType = "PDB", std::string outputFileNaming = "structure"
void Carbohydrate::Generate3DStructureFiles(std::string fileOutputDirectory, std::string outputFileNaming)
{ // ToDo exception handling in centralized function for writing pdb/off
    // Build the filename and path, add appropriate suffix later
    try
    {
        std::string PathAndFileName;
        if (fileOutputDirectory == "unspecified") // "unspecified" is the default
        {
            PathAndFileName += "./" + outputFileNaming;
        }
        else
        {
            PathAndFileName += fileOutputDirectory + "/" + outputFileNaming;
        }
        // Pdb file
        std::string completeFileName = PathAndFileName + ".pdb";
        std::ofstream outFileStream;
        outFileStream.open(completeFileName.c_str());
        this->WritePdb(outFileStream);
        outFileStream.close();
        // Off file
        completeFileName = PathAndFileName + ".off";
        outFileStream.open(completeFileName.c_str());
        cds::WriteMoleculeToOffFile(this->getResidues(), outFileStream, this->getName());
        outFileStream.close();
    }
    catch(const std::string &exceptionMessage)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "carbohydrateBuilder class caught this exception message: " + exceptionMessage);
        this->SetStatus("ERROR", exceptionMessage);
    }
    catch (const std::runtime_error &error)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
        this->SetStatus("ERROR", error.what());
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "carbohydrateBuilder class caught a throw that was not anticipated. Curious. Death cometh?");
        this->SetStatus("ERROR", "carbohydrateBuilder caught a throw type that was not anticipated. Pretty please report how you got to this to glycam@gmail.com.");
    }
}

//////////////////////////////////////////////////////////
//                  PRIVATE FUNCTIONS                   //
//////////////////////////////////////////////////////////
void Carbohydrate::ApplyDeoxy(ParsedResidue* deoxyResidue)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Dealing with deoxy for " + deoxyResidue->getName());
    ParsedResidue* residueToBeDeoxified = deoxyResidue->GetParent();
    residueToBeDeoxified->MakeDeoxy(deoxyResidue->GetLink());
    this->deleteResidue(deoxyResidue); // Remove the deoxy derivative now.
    return;
}

void Carbohydrate::MoveAtomsFromPrepResidueToParsedResidue(prep::PrepFile& prepResidues, ParsedResidue* parsedResidue)
{
    cds::Residue* prepResidue = prepResidues.getResidue(this->GetGlycamResidueName(parsedResidue));
    if (prepResidue == nullptr)
    {
        std::string message = "Did not find prep entry for " + parsedResidue->getName() + " with glycam residue code: " + this->GetGlycamResidueName(parsedResidue);
        gmml::log(__LINE__,__FILE__,gmml::ERR, message);
        throw(std::runtime_error(message));
    }
    parsedResidue->setName(prepResidue->getName()); // Need parsedResidue to be called e.g. 0MA and not DManpa1-4. I can see this being an issue now. Perhaps need a "GlycamName" variable?
    parsedResidue->setAtoms(prepResidue->extractAtoms()); // This moves the atoms, i.e. for prepResidue "Moved from objects are left in a valid but unspecified state"
    prepResidues.deleteResidue(prepResidue); // Death to the prepResidue, if there are repeats with the same name, the next search would find the one without atoms.
    //std::cout << "Finished moving atoms from prepResidue to parsed Residue. Adventure awaits! Huzzah!" << std::endl;
    return;
}

void Carbohydrate::DerivativeChargeAdjustment(ParsedResidue* parsedResidue)
{
    gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeChargeAdjustmentLookupContainer lookup;
    std::string adjustAtomName = lookup.GetAdjustmentAtom(parsedResidue->getName());
    adjustAtomName += parsedResidue->GetLinkageName().substr(0,1);
    cds::Atom* atomToAdjust = parsedResidue->GetParent()->FindAtom(adjustAtomName);
    atomToAdjust->setCharge(atomToAdjust->getCharge() + lookup.GetAdjustmentCharge(parsedResidue->getName()));
    // Log it:
    std::stringstream ss;
    ss << "Charge Adjustment.\n" << "    Derivative is " << parsedResidue->getName() << ". Adjusting charge on " << atomToAdjust->getName() << "\n";
    ss << "    Adjusting by: " << lookup.GetAdjustmentCharge(parsedResidue->getName()) << "\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
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

void Carbohydrate::ConnectAndSetGeometry(cds::Residue* childResidue, cds::Residue* parentResidue)
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
        parentAtom = cdsSelections::getNonCarbonHeavyAtomNumbered(parentResidue->getAtoms(), linkageLabel.substr(linkPosition));
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
    GeometryTopology::FindAtomsToMoveSetDistance(parentAtom, childAtom);
    //   Now bond the atoms. This could also set distance?, and angle? if passed to function?
    childAtom->addBond(parentAtom); // parentAtom also connected to childAtom. Fancy.
    logss << "Bonded " << parentResidue->getName() << "@" << parentAtom->getName() << " to " << childResidue->getName() << "@" << childAtomName << std::endl;
    // Angle
    logss << "Setting angles.\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    for (auto &parentAtomNeighbor : parentAtom->getNeighbors())
    {
        if ( (parentAtomNeighbor->getName().at(0) != 'H') && (parentAtomNeighbor != childAtom ) )
        {
            std::cout << "Setting angle between\nparentNeighbor " << parentAtomNeighbor->getName() << " " << parentAtomNeighbor->getCoordinate()->ToString() << "\nparent " << parentAtom->getName() << " " << parentAtom->getCoordinate()->ToString() << "\nand child " << childAtom->getName() << " " << childAtom->getCoordinate()->ToString() << "\nchild residue " << childResidue->getName() << " will move\n";
            GeometryTopology::SetAngle(parentAtomNeighbor->getCoordinate(), parentAtom->getCoordinate(), childAtom->getCoordinate(), constants::DEFAULT_ANGLE, childResidue->getCoordinates());
            break;
        }
    }
    // Yo if you do this here, then the atoms that move in RotatableDihedral class won't include atoms that get added later. You need to trigger an update of that if you want to wiggle later.
    std::cout << "Finding rotatable dihedrals and applying metadata." << std::endl;
    cds::ResidueLinkage& linkage = glycosidicLinkages_.emplace_back(childResidue, parentResidue);
    std::cout << "Setting default shape" << std::endl;
    linkage.SetDefaultShapeUsingMetadata();
    std::cout << "Wiggling what we have" << std::endl;
    //    std::vector<cds::Atom*> allAtomsInCarb = this->getAtoms();
    //    linkage.SimpleWiggleCurrentRotamers(allAtomsInCarb, allAtomsInCarb, 5);
    std::vector<cds::Atom*> childAtoms = childResidue->getAtoms(); // keeps them alive in memory
    std::vector<cds::Atom*> parentAtoms = parentResidue->getAtoms(); // keeps them alive in memory
    linkage.SimpleWiggleCurrentRotamers(childAtoms, parentAtoms, 5);
    std::cout << "Overlaps resolved greedily" << std::endl;
    return;
}

// Gonna choke on cyclic glycans. Add a check for IsVisited when that is required.
void Carbohydrate::FigureOutResidueLinkages(cds::Residue* from_this_residue1, cds::Residue* to_this_residue2)
{
    //MolecularModeling::ResidueVector neighbors = to_this_residue2->GetNode()->GetResidueNeighbors();

    // Additional code to sort neighbors by lowest index.
    // Only required so that numbers match those assigned in condensed sequence class
    // Should not be done this way, need a generic graph structure and then to centralize everything.
    //MolecularModeling::ResidueVector neighbors = selection::SortResidueNeighborsByAcendingConnectionAtomNumber(to_this_residue2->GetNode()->GetResidueNodeConnectingAtoms());
    // End addtional sorting code.
    /* Breath first code */
    // for(auto &neighbor : neighbors)
    // {
    //     if(neighbor->GetIndex() != from_this_residue1->GetIndex()) // If not the previous residue
    //     {
    //         residue_linkages->emplace_back(neighbor, to_this_residue2);
    //     }
    // }
    /* End Breath first code */
    for(auto &neighbor : to_this_residue2->getChildren())
    {
        if(neighbor->getIndex() != from_this_residue1->getIndex())
        {
            glycosidicLinkages_.emplace_back(neighbor, to_this_residue2); // Depth first. For Breath first remove this line, and comment out above.
            this->FigureOutResidueLinkages(to_this_residue2, neighbor);
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

void Carbohydrate::SetDefaultShapeUsingMetadata()
{
    for(auto &linkage : glycosidicLinkages_)
    {
        linkage.SetDefaultShapeUsingMetadata();
    }
    return;
}

void Carbohydrate::ResolveOverlaps()
{
    for(auto &linkage : glycosidicLinkages_)
    {
        std::vector<cds::Atom*> allAtomsInCarb = this->getAtoms();
        linkage.SimpleWiggleCurrentRotamers(allAtomsInCarb, allAtomsInCarb, 5);
    }
    return;
}
