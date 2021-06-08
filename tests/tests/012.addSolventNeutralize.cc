#include "../../gmml/includes/gmml.hpp"
#include "../../gmml/includes/MolecularModeling/assembly.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../gmml/includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "../../gmml/includes/ParameterSet/OffFileSpace/offfileprocessingexception.hpp"
#include <iostream>
#include <string>

int main()
{
//prep
	// Note: Oliver was checking the functionality. It does not yet work as required, but it took a while to figure out how to run the code
    // So this test is just a snapshot of how it's currently working and how I managed to get output.
    // Create the carbohydrate assembly
    PrepFileSpace::PrepFile* prepA = new PrepFileSpace::PrepFile("../dat/prep/GLYCAM_06j-1.prep");
    std::string condensed_sequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    MolecularModeling::Assembly carbAssembly;
    carbAssembly.BuildAssemblyFromCondensedSequence (condensed_sequence, prepA);
    // Add the ions
    std::string ion_library_file = "../dat/CurrentParams/other/atomic_ions.lib"; // Guessed.
    std::string ion_parameter_file = "../dat/frcmod/frcmod.ionsff99_tip3p"; // Guessed by looking in Assembly::AddIon
    // Attempt to Neutralize by adding both. Only one type will be actually added.
    carbAssembly.AddIon("Na+", ion_library_file, ion_parameter_file);
    carbAssembly.AddIon("Cl-", ion_library_file, ion_parameter_file);
    // Create the solvent assembly
    MolecularModeling::Assembly solventAssembly("../dat/lib/tip3pbox.off", gmml::LIB); // Guessed.
    // Add the solvent around the carbohydrate
    double boxExtensionFromSolute = 10.0, closenessBetweenWaters = 1.0;
    std::string solvent_library_file = "../dat/CurrentParams/other/solvents.lib"; // Guessed.
    carbAssembly.AddSolvent(boxExtensionFromSolute, closenessBetweenWaters, &solventAssembly, solvent_library_file);
    // Write out the pdb file.
    PdbFileSpace::PdbFile *outputPdbFile = carbAssembly.BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("012.addSolventNeutralize.pdb");
    std::cout << "This is the end, my only friend, the end." << std::endl;
}
//prep file

