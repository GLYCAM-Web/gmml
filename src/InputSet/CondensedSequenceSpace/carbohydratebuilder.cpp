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
    pt::ptree child1, child2, childchild3, childchildchild4;
    pt::ptree glycosidicLinkageNode;
    pt::ptree rotatableDihedralsNode;
    childchildchild4.put("index", "1");
    childchild3.put("sequence", this->GetSequenceString());
    childchild3.add_child("glycosidicLinkages", childchildchild4);
    child1.put("type", "Sequence");
    child2.add_child("Evaluate", childchild3);
    root.add_child("entity", child1);
    root.add_child("responses", child2);
//    root.put("sequence", this->GetSequenceString());
//    for (auto &linkage : *(this->GetGlycosidicLinkages())) // I get back a pointer to the ResidueLinkageVector so I *() it to the first element
//    {
//        std::stringstream ss;
//        ss << "LinkageIndex" << linkage.GetIndex();
//        // We only care about those linkages with multiple rotamers
//        RotatableDihedralVector rotatableDihedrals = linkage.GetRotatableDihedralsWithMultipleRotamers();
//        for (auto &rotatableDihedral : rotatableDihedrals)
//        {
//            glycosidicLinkageNode.put(ss.str(), "All likely rotamers"); // as I only want it to be there if there are multiple rotamers in the linkage
//            for (auto &metadata : rotatableDihedral.GetMetadata())
//            {
//                std::stringstream ssmeta;
//                ssmeta << "Dihedral" << metadata.number_of_bonds_from_anomeric_carbon_ << "_" << metadata.index_;
//                rotatableDihedralsNode.put(metadata.dihedral_angle_name_, metadata.rotamer_name_);
//                glycosidicLinkageNode.add_child(ssmeta.str(), rotatableDihedralsNode);
//            }
//        }
//    }
//    root.add_child("Glycosidic Linkages:", glycosidicLinkageNode);
    pt::write_json(std::cout, root);
}

        /*
         * struct DihedralAngleData
{
    std::string linking_atom1_ ;
    std::string linking_atom2_ ;
    std::string dihedral_angle_name_ ;
    double default_angle_value_ ;
    double lower_deviation_ ;
    double upper_deviation_ ;
    double weight_;
    std::string rotamer_type_ ; // permutation or conformer
    std::string rotamer_name_ ;
    int number_of_bonds_from_anomeric_carbon_;
    int index_ ; // Used to indicate whether multiple entries are meant to overwrite each other or generate an additional angle
    StringVector residue1_conditions_ ;
    StringVector residue2_conditions_ ;
    std::string atom1_ ;
    std::string atom2_ ;
    std::string atom3_ ;
    std::string atom4_ ;
} ; */

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
