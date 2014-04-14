
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorUnrecognizedResidue::PdbPreprocessorUnrecognizedResidue() {}

PdbPreprocessorUnrecognizedResidue::PdbPreprocessorUnrecognizedResidue(string residue_name, char chain_id, int sequence_number) :
    residue_name_(residue_name), residue_chain_id_(chain_id), residue_sequence_number_(sequence_number) {}

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

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorUnrecognizedResidue::Print(ostream &out)
{
    cout << "Residue name: " << residue_name_
         << ", Chain id: " << residue_chain_id_
         << ", Sequence number: " << residue_sequence_number_
         << endl;
}








