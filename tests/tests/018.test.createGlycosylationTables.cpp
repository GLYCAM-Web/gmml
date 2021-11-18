#include <iostream>
#include <string>
#include <vector>
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/InternalPrograms/glycosylationSiteFinder.hpp"

// To become the test
int main(int argc, char* argv[])
{
    if ( (argc != 2) && (argc != 3) )
    {
        std::cout << "Usage: createGlycosylationTables.exe inputFile.pdb [outputFileName]\n";
        std::cout << "Example: pdb2glycam 1RVX.pdb \n";
        std::exit(EXIT_FAILURE);
    }
    MolecularModeling::Assembly assembly (argv[1], gmml::InputFileType::PDB);
    assembly.BuildStructureByDistance(10); // number of threads to use.
    assembly.GenerateResidueNodesInAssembly();
    glycoproteinBuilder::GlycosylationSiteFinder siteFinder(assembly);
    std::cout << siteFinder.PrintTable();
    std::vector<GlycosylationSiteInfo> tableInfo = siteFinder.GetTable();
//    for (auto &tableElement : tableInfo)
//    {
//        std::cout << tableElement.Print() << "\n";
//    }
    return 0;
}
