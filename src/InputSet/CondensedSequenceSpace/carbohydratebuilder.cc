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

std::string carbohydrateBuilder::GetOfficialSequenceString()
{
    return officialSequenceString_;
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

bool carbohydrateBuilder::GetSequenceIsValid()
{
    return sequenceIsValid_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
//                      FUNCTIONS                       //
//////////////////////////////////////////////////////////

void carbohydrateBuilder::GenerateSingle3DStructure(std::string outputFileNaming)
{
    this->SetDefaultShapeUsingMetadata();
    this->ResolveOverlaps();
    this->Write3DStructureFile("PDB", outputFileNaming);
    return;
}

// Not yet functional. The json is read in, but nothing happens with it.
void carbohydrateBuilder::GenerateRotamers(std::string jsonSelection)
{
    if (!jsonSelection.empty())
        this->ReadUserSelectionsJSON("");
    else
        //Generate all rotamers
    return;
}

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
    j_responses["Evaluate"]["officialSequence"] = this->GetOfficialSequenceString();
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
void carbohydrateBuilder::Write3DStructureFile(std::string type, std::string filename)
{
    // Use type to figure out which type to write, eg. PDB OFFFILE etc.
    PdbFileSpace::PdbFile *outputPdbFile = this->GetAssembly()->BuildPdbFileStructureFromAssembly();
    std::string emilyThePdbFilename = filename + ".pdb";
//    filename += ".pdb";
    outputPdbFile->Write(emilyThePdbFilename);
    // Le hack:
    //void OffFileSpace::OffFile::Write(std::string file_name, int CoordinateIndex, MolecularModeling::Assembly* assembly)
    OffFileSpace::OffFile frankTheOffFile;
    int CoordinateIndex = 0;
    frankTheOffFile.Write(filename + ".off", CoordinateIndex, this->GetAssembly());
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

void carbohydrateBuilder::SetOfficialSequenceString(std::string sequence)
{
    officialSequenceString_ = sequence;
}

void carbohydrateBuilder::SetSequenceIsValid(bool isValid)
{
    sequenceIsValid_ = isValid;
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
    assembly_.SetName("CONDENSEDSEQUENCE");
    this->SetInputSequenceString(inputSequenceString);
   // CondensedSequence condensedSeqence(inputSequenceString); // This is all weird. condensedSequence should be merged into this class eventually
   // if (condensedSeqence.ParseSequenceAndCheckSanity(inputSequenceString))
 //   { // if a valid sequence
        this->SetSequenceIsValid(true);
     //   this->SetOfficialSequenceString(condensedSeqence.BuildLabeledCondensedSequence(CondensedSequence::Reordering_Approach::LONGEST_CHAIN, false));
        PrepFileSpace::PrepFile* prepFile = new PrepFileSpace::PrepFile(inputPrepFilePath);
        assembly_.BuildAssemblyFromCondensedSequence(inputSequenceString, prepFile);
        // So in the above BuildAssemblyFromCondensedSequence code, linkages are generated that are inaccessible to me.
        // Condensed sequence should be separated so I can handle everything here, but it's a mess.
        // In the mean time I must also create linkages at this level, but I can't figure out how to reset linkage IDs.
        // Maybe I can check if both passed in residues are the same, and reset if so? That only happens here I think.
        this->FigureOutResidueLinkagesInGlycan(assembly_.GetResidues().at(0), assembly_.GetResidues().at(0), &glycosidicLinkages_);
 //   }
//    else
//    { // if not a valid sequence
//        //hmmm not sure how to handle that. Need to ask Dan
//        this->SetSequenceIsValid(false);
//    }
    return;
}
