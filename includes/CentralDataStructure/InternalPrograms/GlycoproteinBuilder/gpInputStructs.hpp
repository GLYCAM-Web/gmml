#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP

#include <string>
#include <vector>
#include <iostream>

struct GlycositeInput
{
	// Constructor
    GlycositeInput(std::string proteinResidueId, std::string glycanInputType, std::string glycan) : proteinResidueId_ (proteinResidueId), glycanInputType_ (glycanInputType), glycanInput_ (glycan) {}
	std::string proteinResidueId_ = "";  // E.g. ?_20_A if no chain ID and residue number is 20 and insertion code is A. C_20_? if chain id is C and there is no insertion code.
	std::string glycanInputType_ = "";  	// "Library" if pre-build as a pdb file or "Sequence" if glycam condensed nomenclature
	std::string glycanInput_ = "";		// E.g. Man9 if "Library" glycanInputType. E.g. DGlcpNAcb1-4DGlcpNAcb1-OH if "Sequence".
};

struct GlycoproteinBuilderInputs
{
	std::string workingDirectory_ = "Default";				// Default is to figure out current directory.
	std::string prepFileLocation_ = "Default";				// Default is to figure out install directory + ../dat/prep/GLYCAM_06j-1_GAGS.prep.
	std::string substrateFileName_ = "Undefined";			// Program will throw if left as "Undefined".
	std::string number3DStructures_ = "1";					//ToDo Implement this.
	std::string maxThreads_ = "1"; 							//ToDo Implement this.
	std::string persistCycles_ = "5";
	std::string overlapTolerance_ = "0.1";
	std::string isDeterministic_ = "false";
	std::vector<GlycositeInput> glycositesInputVector_;	// No default, program will throw if uninitialized.
};

namespace GPInputs
{ // ToDo make these structs part of a namespace, then simplify names within namespace.
GlycoproteinBuilderInputs readGPInputFile(std::string workingDirectory, std::string inputFileName);
}

#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP
