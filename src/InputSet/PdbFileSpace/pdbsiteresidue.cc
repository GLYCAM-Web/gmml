#include "../../../includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

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
    Trim(residue_name_);
    if(section.substr(4,1) == " ")
        residue_chain_id_ = ' ';
    else
        residue_chain_id_ = ConvertString<char>(section.substr(4,1));
    if(section.substr(5, 4) == "    ")
        residue_sequence_number_ = iNotSet;
    else
        residue_sequence_number_ = ConvertString<int>(section.substr(5, 4));
    if(section.substr(9,1) == " ")
        residue_insertion_code_ = ' ';
    else
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

int PdbSiteResidue::GetResidueSequenceNumber(){
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
void PdbSiteResidue::Print(ostream &out)
{
    out << "Residue Name: " << residue_name_
        << ", Residue Chain ID: " << residue_chain_id_
        << ", Residue Sequence Number: ";
    if(residue_sequence_number_ != iNotSet)
        out << residue_sequence_number_;
    else
        out << " ";
    out << ", Residue Insertion Code: " << residue_insertion_code_
        << endl << endl;
}
