#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp"

using PdbPreprocessorSpace::PdbPreprocessorHistidineMapping;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorHistidineMapping::PdbPreprocessorHistidineMapping() {}

PdbPreprocessorHistidineMapping::PdbPreprocessorHistidineMapping(
        char residue_chain_id, int residue_sequence_number,
        gmml::PdbPreprocessorHISMapping selected_mapping,
        char residue_insertion_code, char residue_alternate_location) :
    residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number),
    selected_mapping_(selected_mapping), residue_insertion_code_(residue_insertion_code),
    residue_alternate_location_(residue_alternate_location) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorHistidineMapping::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorHistidineMapping::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
gmml::PdbPreprocessorHISMapping PdbPreprocessorHistidineMapping::GetSelectedMapping()
{
    return selected_mapping_;
}
char PdbPreprocessorHistidineMapping::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
char PdbPreprocessorHistidineMapping::GetResidueAlternateLocation()
{
    return residue_alternate_location_;
}
std::string PdbPreprocessorHistidineMapping::GetStringFormatOfSelectedMapping()
{
    switch(selected_mapping_)
    {
        case 1:
            return "HIE";
        case 2:
            return "HIP";
        case 3:
            return "HID";
        default:
            return "";
    }
}
std::string PdbPreprocessorHistidineMapping::GetStringFormatOfMapping(gmml::PdbPreprocessorHISMapping his_mapping)
{
    switch(his_mapping)
    {
        case 1:
            return "HIE";
        case 2:
            return "HIP";
        case 3:
            return "HID";
        default:
            return "";
    }
}
std::vector<std::string> PdbPreprocessorHistidineMapping::GetAllHISMappingAsString()
{
    std::vector<std::string> all_his_mapping_as_string;
    for(int his_mapping = gmml::HIE; his_mapping != gmml::HID; his_mapping++)
        all_his_mapping_as_string.push_back(GetStringFormatOfMapping((gmml::PdbPreprocessorHISMapping)his_mapping));
    return all_his_mapping_as_string;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorHistidineMapping::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorHistidineMapping::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorHistidineMapping::SetSelectedMapping(gmml::PdbPreprocessorHISMapping selected_mapping)
{
    selected_mapping_ = selected_mapping;
}
void PdbPreprocessorHistidineMapping::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}
void PdbPreprocessorHistidineMapping::SetResidueAlternateLocation(char residue_alternate_location)
{
    residue_alternate_location_ = residue_alternate_location;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorHistidineMapping::Print(std::ostream &out)
{
    out << "Chain id: " << residue_chain_id_
         << ", Sequence_number: " << residue_sequence_number_
         << ", Selected mapping: " << GetStringFormatOfSelectedMapping()
         << ", Insertion code: " << GetResidueInsertionCode()
         << ", Alternate location: " << GetResidueAlternateLocation()
         << std::endl;
}
