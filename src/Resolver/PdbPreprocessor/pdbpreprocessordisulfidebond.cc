
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
int PdbPreprocessorDisulfideBond::GetResidueNumber1()
{
    return residue_number_1_;
}
int PdbPreprocessorDisulfideBond::GetResidueNumber2()
{
    return residue_number_2_;
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
void PdbPreprocessorDisulfideBond::SetResidueNumber1(int residue_number_1)
{
    residue_number_1_ = residue_number_1;
}
void PdbPreprocessorDisulfideBond::SetResidueNumber2(int residue_number_2)
{
    residue_number_2_ = residue_number_2;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorDisulfideBond::Print(ostream &out)
{
}






