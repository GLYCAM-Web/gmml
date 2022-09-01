#include "includes/ParameterSet/PrepFile/prepFile.hpp"
#include <fstream>

int main()
{
    std::string prepFilePath = "../dat/prep/GLYCAM_06j-1_GAGS.prep";
    //std::string condensed_sequence = "LIdopAa1-4DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-6DGalpb1-4DGlcpNAc[6S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
//    prep::PrepFile glycamPrepFile(prepFilePath);
//    for ( auto &prepResidue : glycamPrepFile.getResidues() )
//    {
//    	prepResidue->SetConnectivities();
//    }
//    std::cout << "*\n*\n*\n*\n*\n*\n*\n*\n*\n";
    //std::vector<std::string> residuesToLoadFromPrep = {"0GA", "4YB", "4uA", "Cake", "4YA"};
    std::vector<std::string> residuesToLoadFromPrep = {"0GA"};
    prep::PrepFile glycamPrepFileSelect(prepFilePath, residuesToLoadFromPrep);
    for ( auto &prepResidue : glycamPrepFileSelect.getResidues() )
    {
    	prepResidue->SetConnectivities();
    	prepResidue->Generate3dStructure();
    }
    // Need a central place for this:
    std::string outName = "./prepAsPdbFile.pdb";
    std::ofstream outFileStream;
    try
    {
        outFileStream.open(outName.c_str());
    }
    catch(...)
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Output file could not be created:\n" + outName);
        throw std::runtime_error("Output file could not be created:\n" + outName);
    }
    try
    {
    	glycamPrepFileSelect.WritePdb(outFileStream);
    }
    catch(...)
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Error when writing pdbFile class to file:\n" + outName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + outName);
    }
}
