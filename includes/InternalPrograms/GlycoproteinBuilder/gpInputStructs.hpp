#ifndef GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP
#define GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP

#include <string>
#include <vector>
#include <iostream>

struct GlycositeInput
{
	// Constructor
	GlycositeInput(std::string proteinResidueId, std::string glycanInputType, std::string glycan) : proteinResidueId_ (proteinResidueId), glycanInputType_ (glycanInputType), glycanInput_ (glycan) {}
	// Data. No default values
	std::string proteinResidueId_;  // E.g. ?_20 if no chain ID and residue number is 20. C_20 if chain id is C.
	std::string glycanInputType_;  	// "Library" if pre-build as a pdb file or "Sequence" if glycam condensed nomenclature
	std::string glycanInput_;		// E.g. Man9 if "Library" glycanInputType. E.g. DGlcpNAcb1-4DGlcpNAcb1-OH if "Sequence".
};

struct GlycoproteinBuilderInputs // This are all strings to make it easy to convert to from JSON at the gems level.
{
	GlycoproteinBuilderInputs() :  workingDirectory_ ("Default"), prepFileLocation_ ("Default"), substrateFileName_ ("Undefined"),
			number3DStructures_("1"), maxThreads_("1"), persistCycles_ ("5"), overlapTolerance_("0.1"), isDeterministic_("false") {}
	std::string workingDirectory_;						// Default is to figure out current directory.
	std::string prepFileLocation_;						// Default is to figure out install directory + ../dat/prep/GLYCAM_06j-1_GAGS.prep.
	std::string substrateFileName_;						// Default is "Undefined", program will throw if left as "Undefined".
	std::string number3DStructures_;					// Default is 1. //ToDo Implement this.
	std::string maxThreads_; 							// Default is 1. //ToDo Implement this.
	std::string persistCycles_;						    // Default is 5.
	std::string overlapTolerance_;						// Default is 0.1
	std::string isDeterministic_;						// Default is false
	std::vector<GlycositeInput> glycositesInputVector_;	// No default, program will throw if uninitialized.
};

namespace GPInputs
{ // ToDo make these structs part of a namespace, then simplify names within namespace.
GlycoproteinBuilderInputs readGPInputFile(std::string workingDirectory, std::string inputFileName);
}

#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GPINPUTSTRUCTS_HPP
