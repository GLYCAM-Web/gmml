#include "includes/gmml.hpp"
#include <string>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << "tests/inputs/preprocessor_input.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    PdbFileSpace::PdbFile pdbFile(argv[1]);
    PdbPreprocessorSpace::PdbPreprocessor preprocessor(pdbFile);
    preprocessor.ApplyPreprocessingWithTheGivenModelNumber();
    pdbFile.WriteWithTheGivenModelNumber("Processed.pdb");
    return 0;
}
