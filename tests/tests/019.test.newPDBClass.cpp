#include "includes/gmml.hpp"
#include "includes/InputSet/PdbFile/pdbFile.hpp"
#include "includes/Resolver/NewPdbPreprocessor/pdbPreprocessorInputs.hpp"

#include <string>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << "tests/inputs/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    pdb::PdbFile pdbFile(argv[1]);
    std::cout << "pdbFile constructed\n" << std::endl;
    pdb::PreprocessorOptions options;
    pdb::PreprocessorInformation whatever = pdbFile.PreProcess(options);
//    PdbPreprocessorSpace::PdbPreprocessor preprocessor(pdbFile);
//    preprocessor.ApplyPreprocessingWithTheGivenModelNumber();
//    pdbFile.WriteWithTheGivenModelNumber("Processed.pdb");
    return 0;
}
