#include "includes/gmml.hpp"
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "includes/ParameterSet/OffFileSpace/offfileprocessingexception.hpp"
#include <iostream>
#include <string>

int main()
{
//prep
	// Note: Oliver was checking the functionality. It does not yet work as required, but it took a while to figure out how to run the code
    // So this test is just a snapshot of how it's currently working and how I managed to get output.
    // Create the carbohydrate assembly
    PrepFileSpace::PrepFile* prepA = new PrepFileSpace::PrepFile("../dat/prep/GLYCAM_06j-1.prep");
    //std::string condensed_sequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    std::string condensed_sequence = "DNeu5Aca2-6DGalpb1-4DGlcpNAcb1-OH";
    MolecularModeling::Assembly carbAssembly;
    std::cout << "Building carb from assembly" << std::endl;
    carbAssembly.BuildAssemblyFromCondensedSequence (condensed_sequence, prepA);
    // Add the ions
    std::string ion_library_file = "../dat/CurrentParams/other/atomic_ions.lib"; // Guessed.
    std::string ion_parameter_file = "../dat/frcmod/frcmod.ionsff99_tip3p"; // Guessed by looking in Assembly::AddIon
    // Attempt to Neutralize by adding both. Only one type will be actually added.
    std::cout << "Adding ions" << std::endl;
    carbAssembly.AddIon("Na+", ion_library_file, ion_parameter_file);
    carbAssembly.AddIon("Cl-", ion_library_file, ion_parameter_file);
    // Create the solvent assembly
    std::cout << "Creating solvent assembly" << std::endl;
    MolecularModeling::Assembly solventAssembly("../dat/lib/tip3pbox.off", gmml::LIB); // Guessed.
    // Add the solvent around the carbohydrate
    double boxExtensionFromSolute = 10.0, closenessBetweenWaters = 1.0;
    std::string solvent_library_file = "../dat/CurrentParams/other/solvents.lib"; // Guessed.
    std::cout << "Adding solvent" << std::endl;
    carbAssembly.AddSolvent(boxExtensionFromSolute, closenessBetweenWaters, &solventAssembly, solvent_library_file);
    // Write out the pdb file.
    std::cout << "Writing PDB file" << std::endl;
    PdbFileSpace::PdbFile *outputPdbFile = carbAssembly.BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("012.addSolventNeutralize.pdb");
    std::cout << "This is the end, my only friend, the end." << std::endl;
}
//prep file

