#include "../../../includes/InputSet/CondensedSequenceSpace/carbohydratebuilder.h"
#include "../../../includes/GeometryTopology/ResidueLinkages/residue_linkage.h"
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/json_parser.hpp>
#include "../../../includes/External_Libraries/json.hpp"

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

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void carbohydrateBuilder::SetInputSequenceString(std::string sequence)
{
    inputSequenceString_ = sequence;
}

void carbohydrateBuilder::SetOfficialSequenceString(std::string sequence)
{
    officialSequenceString_ = sequence;
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

void carbohydrateBuilder::WriteFile(std::string type, std::string filename)
{
    // Use type to figure out which type to write, eg. PDB OFFFILE etc.
    PdbFileSpace::PdbFile *outputPdbFile = this->GetAssembly()->BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write(filename);
    return;
}


void carbohydrateBuilder::GenerateRotamers()
{

}

// This is the version with the boost library. I gave up on it. Keeping it for a while just in case.
//void carbohydrateBuilder::WriteJSON()
//{
//    namespace pt = boost::property_tree;
//    pt::ptree root;
//    pt::ptree entity, sequence, glycosidicLinkages, responses, dihedral, metadataEntry, likelyRotamer, glycosidicLinkage;
//    std::vector <pt::ptree> ptvMetaEntries, ptvDihedrals, ptvPossibleRotamers, ptvLikelyRotamers;

//    for (auto &linkage : *(this->GetGlycosidicLinkages())) // I get back a pointer to the ResidueLinkageVector so I *() it to the first element
//    {
////        // We only care about those linkages with multiple rotamers
//        RotatableDihedralVector LikelyRotatableDihedrals = linkage.GetRotatableDihedralsWithMultipleRotamers();
//        for (auto &rotatableDihedral : LikelyRotatableDihedrals)
//        {
//            ptvMetaEntries.clear();
//            for (auto &metadata : rotatableDihedral.GetMetadata())
//            {
//                metadataEntry.put("", metadata.rotamer_name_);
//                ptvMetaEntries.push_back(metadataEntry);
//            }
////            std::stringstream dihedralID;
////            dihedralID << "Dihedral" << rotatableDihedral.GetMetadata().at(0).number_of_bonds_from_anomeric_carbon_ << "_" << rotatableDihedral.GetMetadata().at(0).index_;
//            for (auto &entry : ptvMetaEntries)
//            {
//                dihedral.push_back(std::make_pair("", entry));
//                ptvDihedrals.push_back(dihedral);
//            }
//            likelyRotamer.add_child(rotatableDihedral.GetMetadata().at(0).dihedral_angle_name_, dihedral);
//        }
//        if (!LikelyRotatableDihedrals.empty())
//        {
//            glycosidicLinkage.add_child("likelyRotamers", likelyRotamer);
//            std::stringstream linkageID;
//            linkageID << "LinkageIndex" << linkage.GetIndex();
//            glycosidicLinkages.add_child(linkageID.str(), glycosidicLinkage);
//        }
//    }

//    // Construct everything into the JSON thing. I've done my own indendation here to help:
//            sequence.put("sequence", this->GetSequenceString());
//        responses.add_child("Evaluate", sequence);
//        responses.add_child("glycosdicLinkages", glycosidicLinkages);
//        entity.put("type", "Sequence");
//    root.add_child("entity", entity);
//    root.add_child("responses", responses);

//    pt::write_json(std::cout, root);
//}

void carbohydrateBuilder::WriteJSON()
{
    using json = nlohmann::json;
    json root, responses, linkages, entries;
    for (auto &linkage : *(this->GetGlycosidicLinkages())) // I get back a pointer to the ResidueLinkageVector so I *() it to the first element
    {
        std::cout << "linko nameo: " << linkage.GetName() << std::endl;
        RotatableDihedralVector likelyRotatableDihedrals = linkage.GetRotatableDihedralsWithMultipleRotamers();
        for (auto &rotatableDihedral : likelyRotatableDihedrals)
        {
            for (auto &metadata : rotatableDihedral.GetMetadata())
            {
               entries["likelyRotamers"][metadata.dihedral_angle_name_] += (metadata.rotamer_name_);
               entries["possibleRotamers"][metadata.dihedral_angle_name_] += (metadata.rotamer_name_);

            }
        }
        if(!likelyRotatableDihedrals.empty()) // Only want ones with multiple entries. See above call.
        {
            linkages[std::to_string(linkage.GetIndex())] += (entries);
            entries.clear(); // Must do this as some entries match, e.g. likelyRotamers, omg, gt.
        }
    }
    responses["Evaluate"]["glycosidicLinkages"] += (linkages);
    responses["Evaluate"]["inputSequence"] = this->GetInputSequenceString();
    responses["Evaluate"]["officialSequence"] = this->GetOfficialSequenceString();
    root["responses"] += responses;
    root["entity"]["type"] = "sequence";
    std::cout << std::setw(4) << root << std::endl;
    std::cout << "Finito" << std::endl;
    return;
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

//void carbohydrateBuilder::FigureOutResidueLinkageName(Residue_linkage linkage)
//{ // I need this function so I can avoid using condensed sequence class.
//    // This should be linkage.FigureOutName?


//}

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
    this->SetInputSequenceString(condensedSequence);
    this->SetOfficialSequenceString(condensedSequence); // Need to actually do this.
    assembly_.SetName("CONDENSEDSEQUENCE");
    assembly_.BuildAssemblyFromCondensedSequence(condensedSequence, prepA);
    this->FigureOutResidueLinkagesInGlycan(assembly_.GetResidues().at(0), assembly_.GetResidues().at(0), &glycosidicLinkages_);
    this->WriteJSON();

    return;
}
