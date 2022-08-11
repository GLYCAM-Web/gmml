#include "includes/ParameterSet/PrepFile/prepFile.hpp"

int main()
{
    std::string prepFilePath = "../dat/prep/GLYCAM_06j-1_GAGS.prep";
    //std::string condensed_sequence = "LIdopAa1-4DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-6DGalpb1-4DGlcpNAc[6S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    prep::PrepFile glycamPrepFile(prepFilePath);
    for ( auto &prepResidue : glycamPrepFile.getResidues() )
    {
    	prepResidue->SetConnectivities();
    }

}
