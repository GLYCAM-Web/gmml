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
    /* https://github.com/nlohmann/json. See also includes/External_Libraries/json.hpp
     * nlohmann::json is a little funky. If you first declare something like root["Evaluate"] = "string",
     * you can't later do += as it figures out it's an object
     * and not a list. You can declare it to be a list though. json empty_array_explicit = json::array();
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
               j_entries["possibleRotamers"][metadata.dihedral_angle_name_] += (metadata.rotamer_name_);
            }
        }
        if(!likelyRotatableDihedrals.empty()) // Only want ones with multiple entries. See above call.
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
    std::cout << j_root << std::endl;
   // std::cout << std::setw(4) << j_root << std::endl;
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
    this->SetInputSequenceString(condensedSequence);
    this->SetOfficialSequenceString(condensedSequence); // Need to actually do this.
    assembly_.SetName("CONDENSEDSEQUENCE");
    assembly_.BuildAssemblyFromCondensedSequence(condensedSequence, prepA);
    // So in the above BuildAss code, linkages are generated that are inaccessible to me.
    // Condensed sequence should be separated so I can handle everything here, but it's a mess.
    // In the mean time I must also create linkages at this level, so I want to reset the indexes first:
    this->FigureOutResidueLinkagesInGlycan(assembly_.GetResidues().at(0), assembly_.GetResidues().at(0), &glycosidicLinkages_);
    this->WriteJSON();

    return;
}
