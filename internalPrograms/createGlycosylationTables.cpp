#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include <iostream>
#include <string>
#include <vector>

#include "../includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp" // bondAtomsByDistance
#include "../includes/CentralDataStructure/InternalPrograms/glycosylationSiteFinder.hpp"

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: createGlycosylationTables.exe inputFile.pdb\n";
        std::cout << "Example: pdb2glycam 1RVX.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    std::string inputFileName = argv[1];
    pdb::PdbFile inputFile(inputFileName);
    cds::bondAtomsAndResiduesByDistance(inputFile.getResidues());
    glycoproteinBuilder::GlycosylationSiteFinder siteFinder(inputFile.getResidues());
    std::cout << siteFinder.PrintTable();
    // std::vector<GlycosylationSiteInfo> tableInfo = siteFinder.GetTable();
    //    for (auto &tableElement : tableInfo)
    //    {
    //        std::cout << tableElement.Print() << "\n";
    //    }
    return 0;
}
