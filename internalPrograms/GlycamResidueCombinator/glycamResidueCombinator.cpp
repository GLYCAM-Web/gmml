#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/Editors/glycamResidueCombinator.hpp"
#include <fstream>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputPrepFile\n";
        std::cout << "Exmpl: " << argv[0] << "../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep";
        std::exit(1);
    }
    std::string inputFile = argv[1];
    std::cout << "Test 027 Input file is " << inputFile << "\n";
    std::vector<std::string> residuesToLoadFromPrep = {"0GA", "0YB"};
    prep::PrepFile glycamPrepFile(inputFile, residuesToLoadFromPrep);
    for (auto& residue : glycamPrepFile.getResidues())
    {
        std::cout << "Generating those combos boi" << std::endl;
        std::vector<cds::Residue*> newResidues = residueCombinator::generateResidueCombinations(residue);
        newResidues.push_back(residue);
        std::ofstream outFileStream;
        std::string fileName = residue->getName() + ".lib";
        outFileStream.open(fileName.c_str());
        cds::WriteResiduesToOffFile(newResidues, outFileStream);
        outFileStream.close();
    }
}
