
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorHistidineMapping::PdbPreprocessorHistidineMapping() {}

PdbPreprocessorHistidineMapping::PdbPreprocessorHistidineMapping(char residue_chain_id, int residue_sequence_number, PdbPreprocessorHISMapping selected_mapping, char residue_insertion_code) :
    residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number), selected_mapping_(selected_mapping), residue_insertion_code_(residue_insertion_code) {}

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
PdbPreprocessorHISMapping PdbPreprocessorHistidineMapping::GetSelectedMapping()
{
    return selected_mapping_;
}
char PdbPreprocessorHistidineMapping::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
string PdbPreprocessorHistidineMapping::GetStringFormatOfSelectedMapping()
{
    switch(selected_mapping_)
    {
        case 1:
            return "HIE";
        case 2:
            return "HIP";
        case 3:
            return "HID";
    }
}
string PdbPreprocessorHistidineMapping::GetStringFormatOfMapping(PdbPreprocessorHISMapping his_mapping)
{
    switch(his_mapping)
    {
        case 1:
            return "HIE";
        case 2:
            return "HIP";
        case 3:
            return "HID";
    }
}
vector<string> PdbPreprocessorHistidineMapping::GetAllHISMappingAsString()
{
    vector<string> all_his_mapping_as_string;
    for(int his_mapping = HIE; his_mapping != HID; his_mapping++)
        all_his_mapping_as_string.push_back(GetStringFormatOfMapping((PdbPreprocessorHISMapping)his_mapping));
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
void PdbPreprocessorHistidineMapping::SetSelectedMapping(PdbPreprocessorHISMapping selected_mapping)
{
    selected_mapping_ = selected_mapping;
}
void PdbPreprocessorHistidineMapping::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorHistidineMapping::Print(ostream &out)
{
    cout << "Chain id: " << residue_chain_id_
         << ", Sequence_number: " << residue_sequence_number_
         << ", Selected mapping: " << GetStringFormatOfSelectedMapping()
            << ", insertion code: " << GetResidueInsertionCode()
         << cout << endl;
}







