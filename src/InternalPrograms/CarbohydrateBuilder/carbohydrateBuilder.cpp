#include "includes/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp"
#include "includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "includes/CodeUtils/logging.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using CondensedSequenceSpace::carbohydrateBuilder;
using CondensedSequenceSpace::CondensedSequence;

carbohydrateBuilder::carbohydrateBuilder(std::string condensedSequence, std::string prepFilePath)
: assembly_(condensedSequence, prepFilePath)
{
    this->InitializeClass(condensedSequence);
}

//////////////////////////////////////////////////////////
//                       ACCESSORS                      //
//////////////////////////////////////////////////////////

CondensedSequence carbohydrateBuilder::GetCondensedSequence()
{
    return condensedSequence_;
}

std::string carbohydrateBuilder::GetInputSequenceString()
{
    return inputSequenceString_;
}


MolecularModeling::Assembly* carbohydrateBuilder::GetAssembly()
{
    return &assembly_;
}

ResidueLinkageVector* carbohydrateBuilder::GetGlycosidicLinkages()
{
    return &glycosidicLinkages_;
}
//////////////////////////////////////////////////////////
//                      FUNCTIONS                       //
//////////////////////////////////////////////////////////

void carbohydrateBuilder::GenerateSingle3DStructureDefaultFiles(std::string fileOutputDirectory)
{
    this->SetDefaultShapeUsingMetadata();
    this->ResolveOverlaps();
    this->Write3DStructureFile(fileOutputDirectory, "PDB", "structure");
    this->Write3DStructureFile(fileOutputDirectory, "OFFFILE", "structure");
    return;
}
void carbohydrateBuilder::GenerateSingle3DStructureSingleFile(std::string fileOutputDirectory, std::string type, std::string outputFileNaming)
{
    this->SetDefaultShapeUsingMetadata();
    this->ResolveOverlaps();
    this->Write3DStructureFile(fileOutputDirectory, type, outputFileNaming);
    return;
}

void carbohydrateBuilder::GenerateSpecific3DStructure(CondensedSequenceSpace::SingleRotamerInfoVector conformerInfo, std::string fileOutputDirectory)
{
//     std::string linkageIndex; // What Dan is calling linkageLabel. Internal index determined at C++ level and given to frontend to track.
//     std::string linkageName; // Can be whatever the user wants it to be, default to same as index.
//     std::string dihedralName; // omg / phi / psi / chi1 / chi2
//     std::string selectedRotamer; // gg / tg / g- etc
//     std::string numericValue; // user entered 64 degrees. Could be a v2 feature.
    // With a conformer (aka rotamerSet), setting will be different as each rotatable_dihedral will be set to e.g. "A", whereas for linkages
    // with combinatorial rotamers (e,g, phi -g/t, omg gt/gg/tg), we need to set each dihedral as specified, but maybe it will be ok to go 
    // through and find the value for "A" in each rotatable dihedral.. yeah actually it should be fine. Leaving comment for time being.
    try 
    {
        for (auto &rotamerInfo : conformerInfo)
        {
            // std::cout << "linkage: " << rotamerInfo.linkageIndex << " " 
            //     << this->convertIncomingRotamerNamesToStandard(rotamerInfo.dihedralName)
            //     << " being set to " << rotamerInfo.selectedRotamer << std::endl;
            int currentLinkageIndex = std::stoi(rotamerInfo.linkageIndex);
            Residue_linkage *currentLinkage = this->selectLinkageWithIndex(glycosidicLinkages_, currentLinkageIndex);
            std::string standardDihedralName = this->convertIncomingRotamerNamesToStandard(rotamerInfo.dihedralName);
            currentLinkage->SetSpecificShape(standardDihedralName, rotamerInfo.selectedRotamer);
        }
        //this->ResolveOverlaps();
        this->Write3DStructureFile(fileOutputDirectory, "PDB", "structure"); 
        this->Write3DStructureFile(fileOutputDirectory, "OFFFILE", "structure"); 
    } catch (const char* msg) {
        std::cerr << msg << std::endl;
    }
}

int carbohydrateBuilder::GetNumberOfShapes(bool likelyShapesOnly)
{
    int numberOfShapes = 1;
    for(auto &linkage : (*this->GetGlycosidicLinkages()))
    {
        numberOfShapes = (numberOfShapes * linkage.GetNumberOfShapes(likelyShapesOnly));
    }
    return numberOfShapes;
}
// Commenting out for as not being used, and will be confusing later. The front-end calls a differnt function that will build a single, specific rotamer.
// void carbohydrateBuilder::GenerateUpToNRotamers(int maxRotamers)
// {
//     std::cout << "Rotamer Permutator\n";
//     ResidueLinkageVector linkagesOrderedForPermutation = this->SplitLinkagesIntoPermutants(*(this->GetGlycosidicLinkages()));
//     this->generateLinkagePermutationsRecursively(linkagesOrderedForPermutation.begin(), linkagesOrderedForPermutation.end(), maxRotamers);
// }
CondensedSequenceSpace::LinkageOptionsVector carbohydrateBuilder::GenerateUserOptionsDataStruct()
{
    CondensedSequenceSpace::LinkageOptionsVector userOptionsForSequence;
    for (auto &linkage : *(this->GetGlycosidicLinkages())) // I get back a pointer to the ResidueLinkageVector so I *() it to the first element
    {
       // std::cout << "linko nameo: " << linkage.GetName() << std::endl;
        CondensedSequenceSpace::DihedralOptionsVector possibleRotamers, likelyRotamers;
        std::vector<std::string> buffer;
        for (auto &rotatableDihedral : linkage.GetRotatableDihedralsWithMultipleRotamers())
        {
            for (auto &metadata : rotatableDihedral.GetMetadata())
            {
                buffer.push_back(metadata.rotamer_name_);
            }
            possibleRotamers.emplace_back(rotatableDihedral.GetName(), buffer);
            buffer.clear();
            for (auto &metadata : rotatableDihedral.GetLikelyMetadata())
            {
                buffer.push_back(metadata.rotamer_name_); 
            }
            likelyRotamers.emplace_back(rotatableDihedral.GetName(), buffer);
            buffer.clear();
        }   
        // If there are multiple rotamers for this linkage
        if (!linkage.GetRotatableDihedralsWithMultipleRotamers().empty()) 
        {   // Build struct in vector with emplace_back via constructor in struct
            userOptionsForSequence.emplace_back(linkage.GetName(), std::to_string(linkage.GetIndex()),
                                            linkage.GetFromThisResidue1()->GetNumber(),
                                            linkage.GetToThisResidue2()->GetNumber(),
                                            likelyRotamers, possibleRotamers);
        }
    }
    return userOptionsForSequence;
}

std::string carbohydrateBuilder::Print()
{
    std::stringstream logss;
    logss << "CarbohydrateBuilder called using sequence: " << this->GetInputSequenceString() << "\n and contains these Residue_linkages:\n";
    for (auto & linkage : glycosidicLinkages_ )
    {
        logss << linkage.Print();
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return logss.str();
}
//////////////////////////////////////////////////////////
//                   PRIVATE FUNCTIONS                  //
//////////////////////////////////////////////////////////

void carbohydrateBuilder::Write3DStructureFile(std::string fileOutputDirectory, std::string fileType, std::string filename)
{
    // Build the filename and path, add appropriate suffix later
    std::string completeFileName;
    if (fileOutputDirectory != "unspecified") // "unspecified" is the default
    {
        completeFileName += fileOutputDirectory + "/";
    }
    completeFileName += filename;
    // Use type to figure out which type to write, eg. PDB OFFFILE etc.
    if (fileType == "PDB") 
    { // int link_card_direction = -1, int connect_card_existance = 1, int model_index = -1 , bool useInputPDBResidueNumbers = true);
        PdbFileSpace::PdbFile *outputPdbFile = this->GetAssembly()->BuildPdbFileStructureFromAssembly(-1, 0,-1, false);
        completeFileName += ".pdb";
        outputPdbFile->Write(completeFileName);
    }
    if (fileType == "OFFFILE")
    {
    // Le hack:
    //void OffFileSpace::OffFile::Write(std::string file_name, int CoordinateIndex, MolecularModeling::Assembly* assembly)
        OffFileSpace::OffFile frankTheOffFile;
        int CoordinateIndex = 0;
        completeFileName += ".off";
        frankTheOffFile.Write(completeFileName, CoordinateIndex, this->GetAssembly());
    }
    return;
}

void carbohydrateBuilder::SetInputSequenceString(std::string sequence)
{
    inputSequenceString_ = sequence;
}

void carbohydrateBuilder::SetDefaultShapeUsingMetadata()
{
    for(auto &linkage : glycosidicLinkages_)
    {
        linkage.SetDefaultShapeUsingMetadata();
    }
    return;
}


void carbohydrateBuilder::ResolveOverlaps() // Need to consider rotamers.
{
    for(auto &linkage : glycosidicLinkages_)
    {
        linkage.SimpleWiggle(assembly_.GetAllAtomsOfAssembly(), assembly_.GetAllAtomsOfAssembly(), 0.1, 5);
    }
    return;
}

// Gonna choke on cyclic glycans. Add a check for IsVisited when that is required.
void carbohydrateBuilder::FigureOutResidueLinkagesInGlycan(MolecularModeling::Residue *from_this_residue1, MolecularModeling::Residue *to_this_residue2, ResidueLinkageVector *residue_linkages)
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
        if(neighbor->getIndex() != from_this_residue1->GetIndex())
        {
            residue_linkages->emplace_back(neighbor->getDeriviedClass(), to_this_residue2); // Depth first. For Breath first remove this line, and comment out above.
            this->FigureOutResidueLinkagesInGlycan(to_this_residue2, neighbor->getDeriviedClass(), residue_linkages);
        }
    }
    return;
}

void carbohydrateBuilder::InitializeClass(std::string inputSequenceString)
{
    this->SetInputSequenceString(inputSequenceString);
    assembly_.SetName("CONDENSEDSEQUENCE"); // Necessary for off file to load into tleap
    //assembly_.BuildAssemblyFromCondensedSequence(inputSequenceString, &prepFile);
        // So in the above BuildAssemblyFromCondensedSequence code, linkages are generated that are inaccessible to me.
        // Condensed sequence should be separated so I can handle everything here, but it's a mess.
        // In the mean time I must also create linkages at this level, but I can't figure out how to reset linkage IDs.
        // Maybe I can check if both passed in residues are the same, and reset if so? That only happens here I think.
    this->FigureOutResidueLinkagesInGlycan(assembly_.GetResidues().at(0), assembly_.GetResidues().at(0), &glycosidicLinkages_);
    this->resetLinkageIDsToStartFromZero(glycosidicLinkages_); /* just a fudge until I figure out how to have linkage ids be sensible
        When you instantiate a condensedSequence it generates a 3D structure, and sets default torsions using Residue_Linkage. That class is decoupled 
        from this class as it needs to be replaced, but for now I'm using both and Residue_Linkages are created in that class, so when they
        are created again via this class, their index numbers are "too high" as they are static variables.
        */ 
    return;
}

/* This is a straight copy from glycoprotein_builder. I need a high level class that deals with both Residue_linkages, ring shapes
 * etc. That way I can create X shapes of a molecule. For now this will do to figure out some implementation details like file naming.
 */
ResidueLinkageVector carbohydrateBuilder::SplitLinkagesIntoPermutants(ResidueLinkageVector &inputLinkages)
{
    ResidueLinkageVector sortedLinkages;
    for(auto &linkage : inputLinkages)
    {
        if(linkage.CheckIfConformer())
        {
            sortedLinkages.push_back(linkage);
        }
        else // if not a conformer
        {
            std::vector<Rotatable_dihedral> rotatableDihedrals = linkage.GetRotatableDihedralsWithMultipleRotamers(); // only want the rotatabe dihedrals within a linkage that have multiple rotamers. Some bonds won't.
            for(auto &rotatableDihedral : rotatableDihedrals)
            {
                Residue_linkage splitLinkage = linkage; // Copy it to get correct info into class
                std::vector<Rotatable_dihedral> temp = {rotatableDihedral};
                splitLinkage.SetRotatableDihedrals(temp);
                sortedLinkages.push_back(splitLinkage);
                //std::cout << "Split out " << splitLinkage.GetFromThisResidue1()->GetId() << "-" << splitLinkage.GetToThisResidue2()->GetId() << " rotamer with number of shapes: " << rotatableDihedral.GetNumberOfRotamers() << "\n";
            }
        }
    }
    return sortedLinkages;
}

// Adapted from resolve_overlaps. For the website, this is actually handled at the gems level
// Goes deep and then as it falls out of the iteration it's setting and writing the shapes. 
void carbohydrateBuilder::generateLinkagePermutationsRecursively(ResidueLinkageVector::iterator linkage, ResidueLinkageVector::iterator end, int maxRotamers, int rotamerCount)
{
    for(int shapeNumber = 0; shapeNumber < linkage->GetNumberOfShapes(); ++shapeNumber)
    {
        ++rotamerCount;
        if (rotamerCount <= maxRotamers)
        {
            linkage->SetSpecificShapeUsingMetadata(shapeNumber);
        //std::cout << linkage->GetFromThisResidue1()->GetId() << "-" << linkage->GetToThisResidue2()->GetId() << ": " << (shapeNumber + 1) << " of " << linkage->GetNumberOfShapes() <<  "\n";
            if(std::next(linkage) != end)
            {
                this->generateLinkagePermutationsRecursively(std::next(linkage), end, maxRotamers, rotamerCount);
            }
            else // At the end
            {
            //Check for issues? Resolve
            //Figure out name of file: http://128.192.9.183/eln/gwscratch/2020/01/10/succinct-rotamer-set-labeling-for-sequences/
            //Write PDB file
                this->Write3DStructureFile("PDB", std::to_string(rotamerCount));
            }
        }
    }
    return;
}

// This could be elsewhere, like in a selection class or a dealing with residue linkages class.
Residue_linkage* carbohydrateBuilder::selectLinkageWithIndex(ResidueLinkageVector &inputLinkages, int indexQuery)
{
    for(auto &linkage : inputLinkages)
    {
        if (linkage.GetIndex() == indexQuery)
            return &linkage;
    }
    throw "Linkage not found in carbohydrateBuilder::selectLinkageWithIndex()";
}

// Just a placeholder until we have a map for linkage ids so the user won't see these underlying ones.
void carbohydrateBuilder::resetLinkageIDsToStartFromZero(ResidueLinkageVector &inputLinkages)
{
    unsigned long long newIndex = 0;
    for(auto &linkage : inputLinkages)
    {
        linkage.SetIndex(newIndex);
        ++newIndex;
    }
}

// In too much of a rush to do this properly, so I'll make it private and 
// so dumb that you'll have to write a proper one. Yes you!
// Metadata belongs in metadata, not in code.
// Should have put this at the point that the name is compared in 
// bool Rotatable_dihedral::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
// Or somehow that comparison should be moved into metadata. I wonder could you overload the 
// == comparison operator for the metadata struct to check through this list when checking the name? 
std::string carbohydrateBuilder::convertIncomingRotamerNamesToStandard(std::string incomingName)
{   // Lamda function to see if string is in the passed in vector
    auto isAinList = [](std::vector<std::string> names, std::string query)
    {
        if (std::find(names.begin(), names.end(), query) != names.end())
            return true;
        return false;
    }; // Thought about converting incomingName to lowercase first.
    if (isAinList({"Omega", "omega", "Omg", "omg", "OMG", "omga", "omg1", "Omg1", "o"}, incomingName))
        return "Omg";
    if (isAinList({"Phi", "phi", "PHI", "h"}, incomingName))
        return "Phi";
    if (isAinList({"Psi", "psi", "PSI", "s"}, incomingName))
        return "Psi";
    if (isAinList({"Chi1", "chi1", "CHI1", "c1"}, incomingName))
        return "Chi1";
    if (isAinList({"Chi2", "chi2", "CHI2", "c2"}, incomingName))
        return "Chi2";
    if (isAinList({"Omega7", "omega7", "OMG7", "omga7", "Omg7", "omg7", "o7"}, incomingName))
        return "Omg7";
    if (isAinList({"Omega8", "omega8", "OMG8", "omga8", "Omg8", "omg8", "o8"}, incomingName))
        return "Omg8";
    if (isAinList({"Omega9", "omega9", "OMG9", "omga9", "Omg9", "omg9", "o9"}, incomingName))
        return "Omg9";
    std::stringstream ss;
    ss << "Specified rotamer name: \"" << incomingName << "\", is not recognized in convertIncomingRotamerNamesToStandard.";
    throw ss.str();
}
