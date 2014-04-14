
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorDisulfideBond::PdbPreprocessorDisulfideBond() {}

PdbPreprocessorDisulfideBond::PdbPreprocessorDisulfideBond(char residue_chain_id_1, char residue_chain_id_2, int residue_sequence_number_1, int residue_sequence_number_2, double distance, bool is_bonded) :
    residue_chain_id_1_(residue_chain_id_1), residue_chain_id_2_(residue_chain_id_2), residue_sequence_number_1_(residue_sequence_number_1), residue_sequence_number_2_(residue_sequence_number_2), distance_(distance), is_bonded_(is_bonded) {}

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
double PdbPreprocessorDisulfideBond::GetDistance()
{
    return distance_;
}
bool PdbPreprocessorDisulfideBond::GetIsBonded()
{
    return is_bonded_;
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
void PdbPreprocessorDisulfideBond::SetDistance(double distance)
{
    distance_ = distance;
}
void PdbPreprocessorDisulfideBond::SetIsBonded(bool is_bonded)
{
    is_bonded_ = is_bonded;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorDisulfideBond::Print(ostream &out)
{
}






