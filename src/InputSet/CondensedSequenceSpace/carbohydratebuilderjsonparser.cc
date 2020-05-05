#include "../../../includes/InputSet/CondensedSequenceSpace/carbohydratebuilder.hpp"
#include "../../../includes/External_Libraries/json.hpp"

/* This is just for reading in JSON. Perhaps would have been better at the GEMS level, and have the struct filled in there
 * Anyways, it's wonky but it works. I need to use this data structure userSelection to generate rotamers.
 * */

namespace CondensedSequenceSpace
{
struct userSelection
{
    std::string linkageLabel;
    std::string dihedralName;
    std::string dihedralUnits;
    std::string dihedralValue;
};
}

using CondensedSequenceSpace::carbohydrateBuilder;

// Can go down through arrays of arrays and find a key.
nlohmann::json findJsonKeyRecursive(const nlohmann::json& j, std::string query_key)
{
    if(j.is_object())
    {
        for(const auto& el: j.items())
        {
            if(el.key() == query_key)
            {
                return el.value();
            }
            else if(el.value().is_object() || el.value().is_array())
            {
                return findJsonKeyRecursive(el.value(), query_key);
            }
        }
    }
    else if(j.is_array())
    {
        for(const auto& item: j)
        {
            return findJsonKeyRecursive(item, query_key);
        }
    }
    return nullptr;
}

void carbohydrateBuilder::ReadUserSelectionsJSON(std::string jsonInput)
{
    using json = nlohmann::json;

    if (jsonInput.empty()) {
    // Temporary/example until jsonInput is set and this is called from elsewhere.
auto exampleJSON = R"({
    "serviceContext": "oligosaccharideModelingServices",
    "object": "linkageGeometries",
    "sessionId": "lengthyUUIDString",
    "originalSequence": "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH",
    "officialSequence": "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH",
    "jobTokenId": "anotherId",
    "linkageGeometries": [{
            "linkageLabel": "1",
            "dihedralName": "omega",
            "dihedralUnits": "rotamerName",
            "dihedralValue": "gg"
        },
        {
            "linkageLabel": "4",
            "dihedralName": "phi",
            "dihedralValue": "1.27749",
            "dihedralUnits": "radians"
        },
        {
            "linkageLabel": "4",
            "dihedralName": "omega",
            "dihedralValue": "gt",
            "dihedralUnits": "rotamerName"
          }
    ]
})"_json;
    jsonInput = exampleJSON.dump(); // I assume jsonInput will be serialized.
    }
    // REMOVE THE ABOVE IF FINISHED TESTING
    json jsonBourne = json::parse(jsonInput);
    json linkageGeometries = findJsonKeyRecursive(jsonBourne, "linkageGeometries");

    std::vector <json> values;
    CondensedSequenceSpace::userSelection steveie;
    for (auto & el : linkageGeometries.items())
    {
 //       std::cout << el.key() << " k:v " << el.value() << "\n";
        values.emplace_back(el.value());
    }
    for (auto & el : values)
    {
        steveie.dihedralName = el["dihedralName"].get<std::string>();
        steveie.linkageLabel = el["linkageLabel"].get<std::string>();
        steveie.dihedralUnits = el["dihedralUnits"].get<std::string>();
        steveie.dihedralValue = el["dihedralValue"].get<std::string>();
 //     std::cout << "steve: " << steveie.dihedralName << ", " << steveie.linkageLabel << ", " << steveie.dihedralUnits << ", " << steveie.dihedralValue << std::endl;
    }
}


//// Allows implicit converstion
//namespace carbohydrateBuilder
//{
//    using json = nlohmann::json;
////    void to_json(json& j, const linkageSelections& ls)
////    {
////        j = json{{"name", ls.name}, {"address", ls.address}, {"age", ls.age}};
////    }

//    void from_json(const json& j, linkageSelections& ls)
//    {
//        j.at("dihedralName").get_to(ls.dihedralName);
//        j.at("linkageLabel").get_to(ls.linkageLabel);
//        j.at("dihedralUnits").get_to(ls.dihedralUnits);
//        j.at("dihedralValue").get_to(ls.dihedralValue);
//    }
//} // namespace

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
