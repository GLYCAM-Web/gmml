#include "includes/InternalPrograms/io.hpp"
#include "includes/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"

int main(int argc, char* argv[])
{
    std::string workingDirectory = Find_Program_workingDirectory(); // Default behaviour.
    std::string inputFile = "input.txt";
    if (argc == 1 || argc > 3)
    {
    	std::cout << "Usage: " << argv[0] << " [inputFile] " << " [relative working directory]\n";
    	std::cout << "The arguments in [] are optional. Assumptions are made if none are provided.\n";
    	std::cout << "Exmpl: " << argv[0] << " input.txt " << " tests/simple\n";
    }
    if (argc == 2)
    {
    	inputFile = argv [1];
    }
    if (argc == 3)
    {
    	inputFile = argv [1];
        workingDirectory = argv[2];
    }
    std::cout << "Working directory is " << workingDirectory << "\n";
    std::cout << "Input file is " << inputFile << "\n";
    GlycoproteinBuilderInputs inputStruct = GPInputs::readGPInputFile(workingDirectory, inputFile);
    GlycoproteinBuilder glycoproteinBuilder(inputStruct);
    glycoproteinBuilder.ResolveOverlaps(); // Default randomize value is true, and output isn't deterministic.
    glycoproteinBuilder.WriteOutputFiles();
    std::cout << "Program got to end ok" << std::endl;
    return 0;
}
