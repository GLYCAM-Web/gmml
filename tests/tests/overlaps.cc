#include "../../includes/gmml.hpp"
#include <string>

using namespace gmml;
using namespace MolecularModeling;

int main(void) {
    Assembly assemblyA("tests/inputs/overlaps_input1.pdb", PDB);
    assemblyA.BuildStructureByDistance();

    Assembly assemblyB("tests/inputs/overlaps_input2.pdb", PDB);
    assemblyB.BuildStructureByDistance();
    
    double overlap = CalculateAtomicOverlaps(assemblyA.GetAllAtomsOfAssembly(), assemblyB.GetAllAtomsOfAssembly());

    std::cout << "Overlap: " << overlap << std::endl;    

    return 0;
}

