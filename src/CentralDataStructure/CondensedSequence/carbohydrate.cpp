#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06ResidueNameGenerator.hpp"
#include "includes/Abstract/absResidue.hpp" // For the Residue::Type
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepAtom.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepResidue.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp" // serializeAtomNumbers
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include <sstream>
#include <cctype>    // isDigit
#include <algorithm> //  std::erase, std::remove

using cdsCondensedSequence::Carbohydrate;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Carbohydrate::Carbohydrate(std::string inputSequence, std::string prepFilePath) : SequenceManipulator {inputSequence}
{
    this->setName("CONDENSEDSEQUENCE");
    this->ReorderSequence();        // So output is consistent regardless of input order e.g. Fuca1-2[Gala1-3]Glca vs
                                    // Gala1-3[Fuca1-2]Glca. Same 3D structure.
    this->SetIndexByConnectivity(); // For reporting residue index numbers to the user
    prep::PrepFile glycamPrepFileSelect(prepFilePath, this->GetGlycamNamesOfResidues());
    for (auto& cdsResidue : this->getResidues())
    { // Move atoms from prep file into parsedResidues.
        if (cdsResidue->GetType() != cds::ResidueType::Deoxy)
        {
            ParsedResidue* parsedResidue = static_cast<ParsedResidue*>(cdsResidue);
            this->MoveAtomsFromPrepResidueToParsedResidue(glycamPrepFileSelect, parsedResidue);
            if (parsedResidue->GetType() == cds::ResidueType::Derivative)
            { // Deal with adjusting charges for derivatives
                this->DerivativeChargeAdjustment(parsedResidue);
            }
        }
    }
    for (auto& cdsResidue : this->getResidues())
    { // Apply any deoxy
        if (cdsResidue->GetType() == cds::ResidueType::Deoxy)
        {
            this->ApplyDeoxy(static_cast<ParsedResidue*>(cdsResidue));
        }
    }
    // Have atom numbers go from 1 to number of atoms. Note this should be after deleting atoms due to deoxy
    cds::serializeNumbers(this->getAtoms());
    // Set 3D structure
    this->DepthFirstSetConnectivityAndGeometry(this->GetTerminal()); // recurve start with terminal
    // Re-numbering is a hack as indices have global scope and two instances give too high numbers.
    unsigned int linkageIndex = 0;
    // Linkages should be Edges to avoid this as they already get renumbered above.
    // gmml::log(__LINE__, __FILE__, gmml::INF, "Determining what moves in linkages of " + inputSequence);
    for (auto& linkage : glycosidicLinkages_) // These will exist on the vector in order of edge connectivity set above.
    { // Greedy first means the atoms-to-move needs to be updated for every linkage:
        linkage.SetIndex(linkageIndex++);
        //        std::stringstream ss;
        //  ss << linkage.GetName() << " linkage index in ctor is " << linkage.GetIndex();
        //       gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
        linkage.DetermineAtomsThatMove();
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Final carbohydrate overlap resolution starting.");
    this->ResolveOverlaps();
    gmml::log(__LINE__, __FILE__, gmml::INF,
              "Final carbohydrate overlap resolution finished. Returning from carbohydrate ctor");
    return;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void Carbohydrate::deleteResidue(cds::Residue* byeBye)
{ // ToDo Have to do this because ResidueLInkages are not Edges!!! Oliver just make them Edges already.
    for (cds::ResidueLinkage& linkage : glycosidicLinkages_)
    {
        if (linkage.GetFromThisResidue1() == byeBye || linkage.GetToThisResidue2() == byeBye)
        {
            this->deleteLinkage(&linkage);
        }
    }
    cds::Molecule::deleteResidue(byeBye);
}

void Carbohydrate::replaceAglycone(cds::Residue* newAglycone)
{
    for (cds::ResidueLinkage& linkage : glycosidicLinkages_)
    {
        if (linkage.GetFromThisResidue1() == this->GetAglycone() || linkage.GetToThisResidue2() == this->GetAglycone())
        { // Old aglycone atoms are connected to reducing residue, are found during creation of rotatable dihedrals.
            this->cds::Molecule::deleteResidue(this->GetAglycone());
            newAglycone->addNeighbor(newAglycone->getName() + "-" + this->GetReducingResidue()->getName(),
                                     this->GetReducingResidue());
            linkage = cds::ResidueLinkage(this->GetReducingResidue(), newAglycone);
            linkage.SetDefaultShapeUsingMetadata();
            return;
        }
    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
// std::string fileOutputDirectory = "unspecified", std::string fileType = "PDB", std::string outputFileNaming =
// "structure"
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
        this->WriteOff(outFileStream);
        outFileStream.close();
    }
    catch (const std::string& exceptionMessage)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "carbohydrate class caught this exception message: " + exceptionMessage);
        throw exceptionMessage;
    }
    catch (const std::runtime_error& error)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
        throw error;
    }
    catch (...)
    {
        std::string message = "carbohydrate class caught a throw type that was not anticipated. Pretty please report "
                              "how you got to this to glycam@gmail.com.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
}

std::string Carbohydrate::GetNumberOfShapes(bool likelyShapesOnly) const
{
    if (this->CountShapes(likelyShapesOnly) > 4294967296)
    {
        return ">2^32";
    }
    return std::to_string(this->CountShapes(likelyShapesOnly));
}

unsigned long int Carbohydrate::CountShapes(bool likelyShapesOnly) const
{
    unsigned long long int numberOfShapes = 1;
    for (auto& linkage : glycosidicLinkages_)
    {
        numberOfShapes = (numberOfShapes * linkage.GetNumberOfShapes(likelyShapesOnly));
    }
    return numberOfShapes;
}

cds::Residue* Carbohydrate::GetReducingResidue()
{ // Kindly dumb, but centralized stupidity. Tagging during construction would be better. Ano-ano linkages won't work,
  // but this is used by gp builder so ok.
    for (auto residue : this->getResidues())
    { // Return the first sugar residue that isn't an Aglycone
        if (residue->GetType() != ResidueType::Aglycone && residue->GetType() == ResidueType::Sugar)
        {
            // gmml::log(__LINE__, __FILE__, gmml::INF, "Assuming that reducing residue is " + residue->getStringId());
            return residue;
        }
    }
    if (this->GetResidueCount() > 1)
    {
        //        gmml::log(__LINE__, __FILE__, gmml::INF,
        //                  "No ResidueType's assigned, so assuming that reducing residue is " +
        //                      this->getResidues().at(1)->getStringId());
        return this->getResidues().at(1);
    }
    std::string message = "Reducing residue requested for Carbohydrate with name " + this->getName() +
                          ", but it doesn't have more than 1 residue";
    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
    throw std::runtime_error(message);
}

cds::Residue* Carbohydrate::GetAglycone()
{ // Kindly dumb, but centralized stupidity. Tagging during construction would be better. Ano-ano linkages won't work,
  // but this is used by gp builder so ok.
    for (auto residue : this->getResidues())
    {
        if (residue->GetType() == ResidueType::Aglycone)
        {
            return residue;
        }
    }
    if (this->GetResidueCount() > 0)
    {

        //        gmml::log(__LINE__, __FILE__, gmml::INF,
        //                  "A ResidueType::Aglycone residue was not present, so assuming that aglycone residue is " +
        //                      this->getResidues().front()->getStringId());
        return this->getResidues().front();
    }
    std::string message = "Aglycone residue requested for Carbohydrate with name " + this->getName() +
                          ", but it doesn't have even 1 residue";
    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
    throw std::runtime_error(message);
}

cds::Atom* Carbohydrate::GetAnomericAtom()
{
    return cdsSelections::guessAnomericAtom(this->GetReducingResidue());
}

//////////////////////////////////////////////////////////
//                  PRIVATE MUTATOR                     //
//////////////////////////////////////////////////////////
void Carbohydrate::deleteLinkage(cds::ResidueLinkage* linkage)
{
    glycosidicLinkages_.erase(std::remove(glycosidicLinkages_.begin(), glycosidicLinkages_.end(), *linkage),
                              glycosidicLinkages_.end());
}

//////////////////////////////////////////////////////////
//                  PRIVATE FUNCTIONS                   //
//////////////////////////////////////////////////////////
void Carbohydrate::ApplyDeoxy(ParsedResidue* deoxyResidue)
{
    // gmml::log(__LINE__, __FILE__, gmml::INF, "Dealing with deoxy for " + deoxyResidue->getName());
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
        std::string message = "Did not find prep entry for " + parsedResidue->getName() +
                              " with glycam residue code: " + this->GetGlycamResidueName(parsedResidue);
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw(std::runtime_error(message));
    }
    parsedResidue->setName(
        prepResidue->getName()); // Need parsedResidue to be called e.g. 0MA and not DManpa1-4. I can see this being an
                                 // issue now. Perhaps need a "GlycamName" variable?
    parsedResidue->setAtoms(prepResidue->extractAtoms()); // This moves the atoms, i.e. for prepResidue "Moved from
                                                          // objects are left in a valid but unspecified state"
    prepResidues.deleteResidue(prepResidue); // Death to the prepResidue, if there are repeats with the same name, the
                                             // next search would find the one without atoms.
    // std::cout << "Finished moving atoms from prepResidue to parsed Residue. Adventure awaits! Huzzah!" << std::endl;
    return;
}

void Carbohydrate::DerivativeChargeAdjustment(ParsedResidue* parsedResidue)
{
    std::string adjustAtomName = GlycamMetadata::GetAdjustmentAtom(parsedResidue->getName());
    adjustAtomName             += parsedResidue->GetLinkageName().substr(0, 1);
    cds::Atom* atomToAdjust    = parsedResidue->GetParent()->FindAtom(adjustAtomName);
    atomToAdjust->setCharge(atomToAdjust->getCharge() + GlycamMetadata::GetAdjustmentCharge(parsedResidue->getName()));
    // Log it:
    //    std::stringstream ss;
    //    ss << "Charge Adjustment.\n"
    //       << "    Derivative is " << parsedResidue->getName() << ". Adjusting charge on " << atomToAdjust->getName()
    //       << "\n";
    //    ss << "    Adjusting by: " << lookup.GetAdjustmentCharge(parsedResidue->getName()) << "\n";
    //    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    return;
}

void Carbohydrate::EnsureIntegralCharge(double charge)
{
    //    std::stringstream ss;
    //    ss << std::fixed;
    //    ss << "Total charge is: " << std::setprecision(5) << charge << std::endl;
    //    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
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
    using cds::Atom;
    using cds::ResidueType;
    std::string linkageLabel = static_cast<ParsedResidue*>(childResidue)->GetLinkageName();
    //    gmml::log(__LINE__, __FILE__, gmml::INF,
    //              "Here with child " + childResidue->getStringId() + " and parent: " + parentResidue->getStringId() +
    //                  " and linkageLabel: " + linkageLabel);
    // This is using the new Node<Residue> functionality and the old AtomNode
    // Now go figure out how which Atoms to bond to each other in the residues.
    // Rule: Can't ever have a child aglycone or a parent derivative.
    Atom* parentAtom         = nullptr;
    std::string childAtomName, parentAtomName;
    if (parentResidue->GetType() == ResidueType::Aglycone)
    {
        parentAtomName = GlycamMetadata::GetConnectionAtomForResidue(parentResidue->getName());
        parentAtom     = parentResidue->FindAtom(parentAtomName);
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
            std::string message =
                "The deduced linkageLabel is too small:\n" + linkageLabel +
                ".\nWe require anomer, start atom number, a dash, and connecting atom number. Example:\na1-4";
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        if (!isdigit(linkageLabel.substr(linkPosition).at(0)))
        {
            std::string message = "Could not convert the last linkage number to an integer: " + linkageLabel;
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        parentAtom =
            cdsSelections::getNonCarbonHeavyAtomNumbered(parentResidue->getAtoms(), linkageLabel.substr(linkPosition));
        // std::cout << "Parent atom is " << parentAtom->getName() << std::endl;
    }
    else
    {
        std::string message = "Error: parent residue: " + parentResidue->getName() +
                              " isn't either Aglycone or Sugar, and derivatives cannot be parents.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    if (parentAtom == nullptr)
    {
        std::string message = "Did not find connection atom in residue: " + parentResidue->getStringId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    // Now get child atom
    if (childResidue->GetType() == ResidueType::Derivative)
    {
        std::string glycamNameForResidue = this->GetGlycamResidueName(static_cast<ParsedResidue*>(childResidue));
        childAtomName                    = GlycamMetadata::GetConnectionAtomForResidue(glycamNameForResidue);
    }
    else if (childResidue->GetType() == ResidueType::Sugar)
    {
        std::string childLinkageNumber = linkageLabel.substr(1, 1);
        if (!isdigit(childLinkageNumber.at(0)))
        {
            std::string message = "Could not convert the first linkage number to an integer: " + childLinkageNumber;
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        childAtomName = "C" + childLinkageNumber;
    }
    else
    {
        std::string message = "Error: child residue: " + childResidue->getName() +
                              " is neither derivative or Sugar (aglycones cannot be children)";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    Atom* childAtom = childResidue->FindAtom(childAtomName);
    if (childAtom == nullptr)
    {
        std::string message =
            "Did not find child atom named " + childAtomName + " in child residue: " + childResidue->getStringId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    // Geometry
    cds::FindAtomsToMoveSetDistance(parentAtom, childAtom);
    //   Now bond the atoms. This could also set distance?, and angle? if passed to function?
    childAtom->addBond(parentAtom); // parentAtom also connected to childAtom. Fancy.
    for (auto& parentAtomNeighbor : parentAtom->getNeighbors())
    {
        if ((parentAtomNeighbor->getName().at(0) != 'H') && (parentAtomNeighbor != childAtom))
        {
            // std::cout << "Setting angle between\nparentNeighbor " << parentAtomNeighbor->getName() << " " <<
            // parentAtomNeighbor->getCoordinate()->ToString() << "\nparent " << parentAtom->getName() << " " <<
            // parentAtom->getCoordinate()->ToString() << "\nand child " << childAtom->getName() << " " <<
            // childAtom->getCoordinate()->ToString() << "\nchild residue " << childResidue->getName() << " will
            // move\n";
            cds::SetAngle(parentAtomNeighbor->getCoordinate(), parentAtom->getCoordinate(), childAtom->getCoordinate(),
                          constants::DEFAULT_ANGLE, childResidue->getCoordinates());
            break;
        }
    }
    // GREEDY: taken care of, but note that the atoms that move in RotatableDihedral class need to be updated after more
    // residues are added.
    cds::ResidueLinkage& linkage = glycosidicLinkages_.emplace_back(childResidue, parentResidue);
    linkage.SetDefaultShapeUsingMetadata();
    std::vector<cds::Atom*> childAtoms  = childResidue->getAtoms();  // keeps them alive in memory
    std::vector<cds::Atom*> parentAtoms = parentResidue->getAtoms(); // keeps them alive in memory
    linkage.SimpleWiggleCurrentRotamers(childAtoms, parentAtoms, 5);
    return;
}

// Gonna choke on cyclic glycans. Add a check for IsVisited when that is required.
void Carbohydrate::DepthFirstSetConnectivityAndGeometry(cds::Residue* currentParent)
{
    // MolecularModeling::ResidueVector neighbors = to_this_residue2->GetNode()->GetResidueNeighbors();

    // Additional code to sort neighbors by lowest index.
    // Only required so that numbers match those assigned in condensed sequence class
    // Should not be done this way, need a generic graph structure and then to centralize everything.
    // MolecularModeling::ResidueVector neighbors =
    // selection::SortResidueNeighborsByAcendingConnectionAtomNumber(to_this_residue2->GetNode()->GetResidueNodeConnectingAtoms());
    // End addtional sorting code.
    // Breath first code
    // for(auto &neighbor : neighbors)
    // {
    //     if(neighbor->GetIndex() != from_this_residue1->GetIndex()) // If not the previous residue
    //     {
    //         residue_linkages->emplace_back(neighbor, to_this_residue2);
    //     }
    // }
    // End Breath first code
    for (auto& child : currentParent->getChildren())
    {
        // glycosidicLinkages_.emplace_back(neighbor, to_this_residue2); // Depth first. For Breath first remove this
        // line, and comment out above.
        //        std::cout << "Setting connection between " << child->getName() << " and its parent " <<
        //        currentParent->getName() << ", the connection has linkageLabel: "
        //         << static_cast<ParsedResidue*>(child)->GetLinkageName() << "\n";
        this->ConnectAndSetGeometry(child, currentParent);
        this->DepthFirstSetConnectivityAndGeometry(child);
    }
    return;
}

std::vector<std::string> Carbohydrate::GetGlycamNamesOfResidues() const
{
    std::vector<std::string> names(this->getResidues().size()); // set size of vec for speed.
                                                                //    std::cout << "Glycam names are: ";
    for (auto& residue : this->getResidues())
    {
        if (residue->GetType() != cds::ResidueType::Deoxy)
        {
            names.push_back(this->GetGlycamResidueName(static_cast<ParsedResidue*>(residue)));
            //            std::cout << names.back() << ", ";
        }
    }
    //    std::cout << "\n";
    return names;
}

std::string Carbohydrate::GetGlycamResidueName(ParsedResidue* residue) const
{
    if (residue->GetType() == cds::ResidueType::Deoxy)
    {
        gmml::log(
            __LINE__, __FILE__, gmml::WAR,
            "Bad idea: We asked for Glycam Residue Name of a deoxy type residue (e.g. the 6D of Glc[6D]) with name: " +
                residue->GetResidueName());
        return "";
    }
    std::string linkages = "";
    if (residue->GetType() == cds::ResidueType::Sugar)
    {
        //        gmml::log(__LINE__, __FILE__, gmml::INF,
        //                  "Checking for glycosidic linkages that connect to " + residue->GetResidueName());
        linkages = residue->GetChildLinkagesForGlycamResidueNaming();
    }
    std::string code = GlycamMetadata::Glycam06ResidueNameGenerator(
        linkages, residue->GetIsomer(), residue->GetResidueName(), residue->GetRingType(),
        residue->GetResidueModifier() + residue->GetRingShape(), residue->GetConfiguration());
    return code;
}

void Carbohydrate::SetDefaultShapeUsingMetadata()
{
    for (auto& linkage : glycosidicLinkages_)
    {
        linkage.SetDefaultShapeUsingMetadata();
    }
    return;
}

void Carbohydrate::ResolveOverlaps()
{
    for (auto& linkage : glycosidicLinkages_)
    {
        std::vector<cds::Residue*> carbResidues = this->getResidues();
        linkage.SimpleWiggleCurrentRotamers(linkage.GetNonReducingOverlapResidues(),
                                            linkage.GetReducingOverlapResidues(), 5);
    }
    return;
}
