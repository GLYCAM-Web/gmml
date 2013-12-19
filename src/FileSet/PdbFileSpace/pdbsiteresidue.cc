#include "../../../includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSiteResidue::PdbSiteResidue() {}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbSiteResidue::GetRecordName(){
    return record_name_;
}

char PdbSiteResidue::GetResidueChainId(){
    return residue_chain_id_;
}

int PdbSiteResidue::GetresidueSequenceNumber(){
    return residue_sequence_number_;
}

char PdbSiteResidue::GetResidueInsertionCode(){
    return residue_insertion_code_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSiteResidue::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbSiteResidue::SetResidueChainId(char residue_chain_id){
    residue_chain_id_ = residue_chain_id;
}

void PdbSiteResidue::SetResidueSequenceNumber(int residue_sequence_number){
    residue_sequence_number_ = residue_sequence_number;
}

void PdbSiteResidue::SetResidueInsertionCode(char residue_insertion_code){
    residue_insertion_code_ = residue_insertion_code;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////






