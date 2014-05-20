
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorUnrecognizedResidue::PdbPreprocessorUnrecognizedResidue() {}

PdbPreprocessorUnrecognizedResidue::PdbPreprocessorUnrecognizedResidue(string residue_name, char chain_id, int sequence_number, char residue_insertion_code) :
    residue_name_(residue_name), residue_chain_id_(chain_id), residue_sequence_number_(sequence_number), residue_insertion_code_(residue_insertion_code) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorUnrecognizedResidue::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorUnrecognizedResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
string PdbPreprocessorUnrecognizedResidue::GetResidueName()
{
    return residue_name_;
}
char PdbPreprocessorUnrecognizedResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorUnrecognizedResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorUnrecognizedResidue::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorUnrecognizedResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedResidue::Print(ostream &out)
{
    out << "Residue name: " << residue_name_
         << ", Chain id: " << residue_chain_id_
         << ", Sequence number: " << residue_sequence_number_
            << ", Insertion code: " << residue_insertion_code_
         << endl;
}








