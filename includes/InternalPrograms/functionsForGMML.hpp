#ifndef GMML_INCLUDES_INTERNALPROGRAMS_FUNCTIONSFORGMML_HPP
#define GMML_INCLUDES_INTERNALPROGRAMS_FUNCTIONSFORGMML_HPP
#include "includes/gmml.hpp"

namespace gmml
{
	void WritePDBFile(MolecularModeling::Assembly &ass, std::string workingDirectory, std::string fileNamePrefix, bool includeOutputFileCount = true);
	void WriteOffFile(MolecularModeling::Assembly &ass, std::string workingDirectory, std::string fileNamePrefix, bool includeOutputFileCount = true);
	bool startsWith(std::string bigString, std::string smallString);
    int CountInternalBonds(Assembly &ass);
}
#endif // GMML_INCLUDES_INTERNALPROGRAMS_FUNCTIONSFORGMML_HPP
