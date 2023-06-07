//#include "includes/CodeUtils/directories.hpp"
#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/inputs.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile\n";
        std::cout << "Exmpl: " << argv[0] << " input.txt\n";
        std::exit(1);
    }
    std::string inputFile = argv[1];
    std::cout << "Input file is " << inputFile << "\n";
    gmmlPrograms::WiggleToSiteInputs inputStruct(inputFile);
    std::cout << "Reading input file complete\n" << std::flush;
    std::cout << inputStruct.Print();
    gmmlPrograms::WiggleToSite wiggler(inputStruct);
    std::cout << "Program got to end ok" << std::endl;
    return 0;
}
