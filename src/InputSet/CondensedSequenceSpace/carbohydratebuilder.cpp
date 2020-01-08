#include "../../../includes/InputSet/CondensedSequenceSpace/carbohydratebuilder.h"
#include "../../../includes/GeometryTopology/ResidueLinkages/residue_linkage.h"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using CondensedSequenceSpace::carbohydrateBuilder;

//carbohydrateBuilder::carbohydrateBuilder()
//{
//    this->InitializeClass();
//}

carbohydrateBuilder::carbohydrateBuilder(std::string selectedBuildType, std::string condensedSequence, std::string prepFilePath)
{
    this->InitializeClass(selectedBuildType, condensedSequence, prepFilePath);
}

//////////////////////////////////////////////////////////
//                      FUNCTIONS                       //
//////////////////////////////////////////////////////////

void carbohydrateBuilder::SetDefaultDihedralAngleGeometryWithMetadata()
{
    for(auto &linkage : glycosidicLinkages_)
    {
        linkage.SetDefaultShapeUsingMetadata();
    }
    return;
}

void carbohydrateBuilder::ResolveOverlaps()
{
    for(auto &linkage : glycosidicLinkages_)
    {
        linkage.SimpleWiggle(assembly_.GetAllAtomsOfAssembly(), assembly_.GetAllAtomsOfAssembly(), 0.1, 5);
    }
    return;
}

//////////////////////////////////////////////////////////
//                   PRIVATE FUNCTIONS                  //
//////////////////////////////////////////////////////////

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

void carbohydrateBuilder::InitializeClass(std::string selectedBuildType, std::string condensedSequence, std::string prepFilePath)
{
    PrepFileSpace::PrepFile* prepA = new PrepFileSpace::PrepFile(prepFilePath);
   // assembly_ = MolecularModeling::Assembly();
    assembly_.SetName("CONDENSEDSEQUENCE");
    assembly_.BuildAssemblyFromCondensedSequence(condensedSequence, prepA);
    this->FigureOutResidueLinkagesInGlycan(assembly_.GetResidues().at(0), assembly_.GetResidues().at(0), &glycosidicLinkages_);
    return;
}
