#include "../../../includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSiteResidue::PdbSiteResidue() {}
PdbSiteResidue::PdbSiteResidue(const string &residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code)
    :residue_name_(residue_name), residue_chain_id_(residue_chain_id), residue_sequence_number_(residue_sequence_number), residue_insertion_code_(residue_insertion_code){}

PdbSiteResidue::PdbSiteResidue(const string& section)
{
    residue_name_ = section.substr(0, 3);
    residue_chain_id_ = ConvertString<char>(section.substr(4,1));
    residue_sequence_number_ = ConvertString<int>(section.substr(5, 4));
    residue_insertion_code_ = ConvertString<char>(section.substr(9,1));
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbSiteResidue::GetResidueName(){
    return residue_name_;
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

void PdbSiteResidue::SetResidueName(const string residue_name){
    residue_name_ = residue_name;
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






