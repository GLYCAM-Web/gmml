#include "includes/CentralDataStructure/Writers/cdsOffWriter.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepPermutator.hpp"
#include "includes/CodeUtils/logging.hpp"
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
    std::vector<std::string> residuesToLoadFromPrep = {"0GA", "4YB", "4uA", "Cake", "4YA"};
    //std::vector<std::string> residuesToLoadFromPrep = {"0GA"};
    prep::PrepFile glycamPrepFileSelect(prepFilePath, residuesToLoadFromPrep);
    // PREP residues
    glycamPrepFileSelect.Write("./prepAsPrepFile.prep");
    // PDB
    std::string fileName = "./prepAsPdbFile.pdb";
    std::ofstream outFileStream;
    try
    {
    	outFileStream.open(fileName.c_str());
    	glycamPrepFileSelect.WritePdb(outFileStream);
    	outFileStream.close();
    }
    catch(...)
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Error when writing pdbFile class to file:\n" + fileName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + fileName);
    }
    // OFF molecule
    glycamPrepFileSelect.setName("MOLECULE");
    fileName = "./prepAsOffFile.off";
    try
    {
    	std::ofstream outFileStream;
    	outFileStream.open(fileName.c_str());
    	cds::WriteMoleculeToOffFile(glycamPrepFileSelect.getResidues(), outFileStream, glycamPrepFileSelect.getName());
    	outFileStream.close();
    }
    catch(...)
    {
    	gmml::log(__LINE__,__FILE__,gmml::ERR, "Error when writing to file:\n" + fileName);
    	throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
    // OFF separate residues
    glycamPrepFileSelect.setName("LIBRARY");
    fileName = "./prepAsLibFile.lib";
    try
    {
    	std::ofstream outFileStream;
    	outFileStream.open(fileName.c_str());
    	cds::WriteResiduesToOffFile(glycamPrepFileSelect.getResidues(), outFileStream);
    	outFileStream.close();
    }
    catch(...)
    {
    	gmml::log(__LINE__,__FILE__,gmml::ERR, "Error when writing to file:\n" + fileName);
    	throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
    // Ok now permutator:
    residuesToLoadFromPrep = {"0GA"};
    prep::PrepFile glycamPreputator(prepFilePath, residuesToLoadFromPrep);
    for(auto & residue : glycamPreputator.getResidues())
    {
        prep::generatePrepPermuations(static_cast<prep::PrepResidue*>(residue));
    }

}
