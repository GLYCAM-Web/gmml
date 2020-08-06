#include "../../../includes/InputSet/CondensedSequenceSpace/carbohydratebuilder.hpp"
#include "../../../includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp"
#include "../../../includes/External_Libraries/json.hpp"
#include "../../../includes/ParameterSet/OffFileSpace/offfile.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using CondensedSequenceSpace::carbohydrateBuilder;
using CondensedSequenceSpace::CondensedSequence;
//carbohydrateBuilder::carbohydrateBuilder()
//{
//    this->InitializeClass();
//}

carbohydrateBuilder::carbohydrateBuilder(std::string condensedSequence, std::string prepFilePath)
{
    this->InitializeClass(condensedSequence, prepFilePath);
}

//////////////////////////////////////////////////////////
//                       ACCESSORS                      //
//////////////////////////////////////////////////////////

CondensedSequence carbohydrateBuilder::GetCondensedSequence()
{
    return condensedSequence_;
}

// std::string carbohydrateBuilder::GetOfficialSequenceString()
// {
//     return officialSequenceString_;
// }

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

// bool carbohydrateBuilder::GetSequenceIsValid()
// {
//     return sequenceIsValid_;
// }

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
//                      FUNCTIONS                       //
//////////////////////////////////////////////////////////

void carbohydrateBuilder::GenerateSingle3DStructure(std::string fileOutputDirectory, std::string type, std::string outputFileNaming)
{
    this->SetDefaultShapeUsingMetadata();
    this->ResolveOverlaps();
    this->Write3DStructureFile(fileOutputDirectory, type, outputFileNaming);
    return;
}

// // Not yet functional. The json is read in, but nothing happens with it.
// void carbohydrateBuilder::GenerateRotamers(std::string jsonSelection)
// {
//     if (!jsonSelection.empty())
//         this->ReadUserSelectionsJSON("");
//     else
//         //Generate all rotamers
//     return;
// }

void carbohydrateBuilder::GenerateRotamer(CondensedSequenceSpace::singleRotamerInfoVector conformerInfo, std::string fileOutputDirectory)
{
//        std::string linkageIndex; // What Dan is calling linkageLabel. Internal index determined at C++ level and given to frontend to track.
//     std::string linkageName; // Can be whatever the user wants it to be, default to same as index.
//     std::string dihedralName; // omg / phi / psi / chi1 / chi2
//     std::string selectedRotamer; // gg / tg / g- etc
//     std::string numericValue; // user entered 64 degrees. Could be a v2 feature.
    // With a conformer (aka rotamerSet), setting will be different as each rotatable_dihedral will be set to e.g. "A", whereas for linkages
    // with combinatorial rotamers (e,g, phi -g/t, omg gt/gg/tg), we need to set each dihedral as specified, but maybe it will be ok to go 
    // through and find the value for "A" in each rotatable dihedral.. yeah actually it should be fine. Leaving comment for time being.
    for (auto &rotamerInfo : conformerInfo)
    {
        int currentLinkageIndex = std::stoi(rotamerInfo.linkageIndex);
        Residue_linkage *currentLinkage = this->selectLinkageWithIndex(glycosidicLinkages_, currentLinkageIndex);
        currentLinkage->SetSpecificShape(rotamerInfo.dihedralName, rotamerInfo.selectedRotamer);
    }
    //this->ResolveOverlaps();
    this->Write3DStructureFile(fileOutputDirectory, "PDB"); // Need to update to pass in output directory
    this->Write3DStructureFile(fileOutputDirectory, "OFFFILE"); // Need to update to pass in output directory

}


void carbohydrateBuilder::GenerateUpToNRotamers(int maxRotamers)
{
    std::cout << "Rotamer Permutator\n";
    ResidueLinkageVector linkagesOrderedForPermutation = this->SplitLinkagesIntoPermutants(*(this->GetGlycosidicLinkages()));
    this->generateLinkagePermutationsRecursively(linkagesOrderedForPermutation.begin(), linkagesOrderedForPermutation.end(), maxRotamers);
}

// This function will be deprecated, JSON to be written at gems level by pydantic
std::string carbohydrateBuilder::GenerateUserOptionsJSON()
{
    /* https://github.com/nlohmann/json. See also includes/External_Libraries/json.hpp
     * nlohmann::json is a little funky. If you first declare something like root["Evaluate"] = "string",
     * you can't later do += as it figures out it's an object
     * and not a list. You can declare it to be a list though. json empty_array_explicit = json::array();
     * Perhaps this would all have been better if I just filled in a struct and then output it as JSON?
     */
    using json = nlohmann::json;
    json j_root, j_responses, j_linkages, j_entries; // using j_ prefix to make clear what is json.
    for (auto &linkage : *(this->GetGlycosidicLinkages())) // I get back a pointer to the ResidueLinkageVector so I *() it to the first element
    {
       // std::cout << "linko nameo: " << linkage.GetName() << std::endl;
        RotatableDihedralVector likelyRotatableDihedrals = linkage.GetRotatableDihedralsWithMultipleRotamers();
        for (auto &rotatableDihedral : likelyRotatableDihedrals)
        {
            for (auto &metadata : rotatableDihedral.GetMetadata())
            {
             j_entries["likelyRotamers"][metadata.dihedral_angle_name_] += (metadata.rotamer_name_);
 // Just a hack for now to get format. I need to call a separate function once I figure out how I'll distinguish All/likely in metadata
             j_entries["possibleRotamers"][metadata.dihedral_angle_name_] += (metadata.rotamer_name_);
         }
     }
        if(!likelyRotatableDihedrals.empty()) // Only want ones with multiple entries. See GetRotatableDihedralsWithMultipleRotamers.
        {   // Order of adding to linkages matters here. I don't know why :(
    j_linkages[std::to_string(linkage.GetIndex())] = (j_entries);
    j_linkages[std::to_string(linkage.GetIndex())]["linkageName"] = linkage.GetName();
            j_entries.clear(); // Must do this as some entries match, e.g. likelyRotamers, omg, gt.
        }
    }
    j_responses["Evaluate"]["glycosidicLinkages"] += (j_linkages);
    j_responses["Evaluate"]["inputSequence"] = this->GetInputSequenceString();
    // j_responses["Evaluate"]["officialSequence"] = this->GetOfficialSequenceString();
    j_responses["Evaluate"]["officialSequence"] = this->GetInputSequenceString();
    j_root["responses"] += j_responses;
    j_root["entity"]["type"] = "sequence";
    // std::cout << j_root << std::endl;
    std::stringstream response;
    response << j_root;
   // std::cout << std::setw(4) << j_root << std::endl;
  //  std::cout << "Finito" << std::endl;
    return response.str();
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
    {
        PdbFileSpace::PdbFile *outputPdbFile = this->GetAssembly()->BuildPdbFileStructureFromAssembly();
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
//void carbohydrateBuilder::Write3DStructureFile(std::string type, std::string filename)
//{
//    // Use type to figure out which type to write, eg. PDB OFFFILE etc.
//    PdbFileSpace::PdbFile *outputPdbFile = this->GetAssembly()->BuildPdbFileStructureFromAssembly();
//    filename += ".pdb";
//    outputPdbFile->Write(filename);
//    return;
//}

void carbohydrateBuilder::SetInputSequenceString(std::string sequence)
{
    inputSequenceString_ = sequence;
}

// void carbohydrateBuilder::SetOfficialSequenceString(std::string sequence)
// {
//     officialSequenceString_ = sequence;
// }

// void carbohydrateBuilder::SetSequenceIsValid(bool isValid)
// {
//     sequenceIsValid_ = isValid;
// }

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

void carbohydrateBuilder::FigureOutResidueLinkagesInGlycan(MolecularModeling::Residue *from_this_residue1, MolecularModeling::Residue *to_this_residue2, ResidueLinkageVector *residue_linkages)
{
    MolecularModeling::ResidueVector neighbors = to_this_residue2->GetNode()->GetResidueNeighbors();
    for(auto &neighbor : neighbors)
    {
        if(neighbor->GetIndex() != from_this_residue1->GetIndex()) // If not the previous residue
        {
            residue_linkages->emplace_back(neighbor, to_this_residue2);
        }
    }
    for(auto &neighbor : neighbors)
    {
        if(neighbor->GetIndex() != from_this_residue1->GetIndex())
        {
            this->FigureOutResidueLinkagesInGlycan(to_this_residue2, neighbor, residue_linkages);
        }
    }
    return;
}

void carbohydrateBuilder::InitializeClass(std::string inputSequenceString, std::string inputPrepFilePath)
{
    // Have to assume that sequence is sane, because using the condensedsequence class functions to check breaks them... probably they should have been private.
    assembly_.SetName("CONDENSEDSEQUENCE");
    this->SetInputSequenceString(inputSequenceString);
    CondensedSequence condensedSeqence(inputSequenceString); // This is all weird. condensedSequence should be merged into this class eventually
    //this->SetOfficialSequenceString(condensedSeqence.BuildLabeledCondensedSequence(CondensedSequence::Reordering_Approach::LONGEST_CHAIN, CondensedSequence::Reordering_Approach::LONGEST_CHAIN, false));
    PrepFileSpace::PrepFile* prepFile = new PrepFileSpace::PrepFile(inputPrepFilePath);
    assembly_.BuildAssemblyFromCondensedSequence(inputSequenceString, prepFile);
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
            RotatableDihedralVector rotatableDihedrals = linkage.GetRotatableDihedralsWithMultipleRotamers(); // only want the rotatabe dihedrals within a linkage that have multiple rotamers. Some bonds won't.
            for(auto &rotatableDihedral : rotatableDihedrals)
            {
                Residue_linkage splitLinkage = linkage; // Copy it to get correct info into class
                RotatableDihedralVector temp = {rotatableDihedral};
                splitLinkage.SetRotatableDihedrals(temp);
                sortedLinkages.push_back(splitLinkage);
                std::cout << "Split out " << splitLinkage.GetFromThisResidue1()->GetId() << "-" << splitLinkage.GetToThisResidue2()->GetId() << " rotamer with number of shapes: " << rotatableDihedral.GetNumberOfRotamers() << "\n";
            }
        }
    }
    return sortedLinkages;
}

//Adapted from resolve_overlaps.
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
    Residue_linkage *linkageFound;
    for(auto &linkage : inputLinkages)
    {
        if (linkage.GetIndex() == indexQuery)
        {
            linkageFound = &linkage;
        }
    }
    return linkageFound;
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