#include "../../../includes/InputSet/CondensedSequenceSpace/carbohydratebuilder.h"
#include "../../../includes/GeometryTopology/ResidueLinkages/residue_linkage.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

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
//                       ACCESSORS                      //
//////////////////////////////////////////////////////////

CondensedSequenceSpace::CondensedSequence carbohydrateBuilder::GetCondensedSequence()
{
    return condensedSequence_;
}

std::string carbohydrateBuilder::GetSequenceString()
{
    return sequenceString_;
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
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void carbohydrateBuilder::SetSequenceString(std::string sequence)
{
    sequenceString_ = sequence;
}

//////////////////////////////////////////////////////////
//                      FUNCTIONS                       //
//////////////////////////////////////////////////////////

void carbohydrateBuilder::GenerateSingle3DStructure()
{
    this->SetDefaultShapeUsingMetadata();
    this->ResolveOverlaps();
    return;
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

void carbohydrateBuilder::GenerateRotamers()
{

}

void carbohydrateBuilder::WriteJSON()
{
    namespace pt = boost::property_tree;
    pt::ptree root;
    root.put("sequence", this->GetSequenceString());
    pt::ptree glycosidicLinkageNode;
    for (auto &linkage : *(this->GetGlycosidicLinkages()))
    {
        std::stringstream li;
        li << "LinkageIndex" << linkage.GetIndex();
        glycosidicLinkageNode.put(li.str(), linkage.GetIndex());
    }
    root.add_child("Glycosidic Linkages:", glycosidicLinkageNode);
    pt::write_json(std::cout, root);
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
    this->SetSequenceString(condensedSequence);
    assembly_.SetName("CONDENSEDSEQUENCE");
    assembly_.BuildAssemblyFromCondensedSequence(condensedSequence, prepA);
    this->FigureOutResidueLinkagesInGlycan(assembly_.GetResidues().at(0), assembly_.GetResidues().at(0), &glycosidicLinkages_);
    this->WriteJSON();

    return;
}
