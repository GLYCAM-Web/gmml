
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorDisulfideBond::PdbPreprocessorDisulfideBond() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorDisulfideBond::GetResidueChainId1()
{
    return residue_chain_id_1_;
}
char PdbPreprocessorDisulfideBond::GetResidueChainId2()
{
    return residue_chain_id_2_;
}
int PdbPreprocessorDisulfideBond::GetResidueSequenceNumber1()
{
    return residue_sequence_number_1_;
}
int PdbPreprocessorDisulfideBond::GetResidueSequenceNumber2()
{
    return residue_sequence_number_2_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorDisulfideBond::SetResidueChainId1(char residue_chain_id_1)
{
    residue_chain_id_1_ = residue_chain_id_1;
}
void PdbPreprocessorDisulfideBond::SetResidueChainId2(char residue_chain_id_2)
{
    residue_chain_id_2_ = residue_chain_id_2;
}
void PdbPreprocessorDisulfideBond::SetResidueSequenceNumber1(int residue_sequence_number_1)
{
    residue_sequence_number_1_ = residue_sequence_number_1;
}
void PdbPreprocessorDisulfideBond::SetResidueSequenceNumber2(int residue_sequence_number_2)
{
    residue_sequence_number_2_ = residue_sequence_number_2;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorDisulfideBond::Print(ostream &out)
{
    cout << "Chain id: " << residue_chain_id_1_
         << ", Sequence number: " << residue_sequence_number_1_
         << ", Chain id: " << residue_chain_id_2_
         << ", Sequence number: " << residue_sequence_number_2_
         << ", Distance: " << distance_
         << endl;
}






