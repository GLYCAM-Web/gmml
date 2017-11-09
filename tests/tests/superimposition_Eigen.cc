#include "../../includes/gmml.hpp"
#include <string>

using namespace gmml;
using namespace MolecularModeling;

int main(void) {
    Assembly moving("tests/inputs/superimposition_Eigen_Moving.pdb", PDB);
    moving.BuildStructureByDistance();

    Assembly target("tests/inputs/superimposition_Eigen_Target.pdb", PDB);
    target.BuildStructureByDistance();
    
    Superimpose(&moving, &target);
    
    PdbFileSpace::PdbFile *outputPdbFile = moving.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile->Write("moved.pdb");

    return 0;
}

