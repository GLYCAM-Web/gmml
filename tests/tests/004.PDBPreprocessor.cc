#include "includes/gmml.hpp"
#include <string>

int main(void)
{
    PdbFileSpace::PdbFile pdbFile("tests/inputs/preprocessor_input.pdb");
    PdbPreprocessorSpace::PdbPreprocessor preprocessor(pdbFile);
    preprocessor.ApplyPreprocessingWithTheGivenModelNumber();
    pdbFile.WriteWithTheGivenModelNumber("Processed.pdb");
    return 0;
}
