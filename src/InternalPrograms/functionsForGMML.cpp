#include "includes/InternalPrograms/functionsForGMML.hpp"

void gmml::WritePDBFile(MolecularModeling::Assembly &ass, std::string workingDirectory, std::string fileNamePrefix, bool includeOutputFileCount)
{
	static int outputFileCount = 0;
	outputFileCount++;
	if (includeOutputFileCount)
	{
		std::stringstream ss;
		ss << "_" << std::setw(5) << std::setfill('0') << outputFileCount;
		fileNamePrefix += ss.str();
	}
	std::string fullOutputFileName = workingDirectory + fileNamePrefix + ".pdb";
	PdbFileSpace::PdbFile *outputPdbFile = ass.BuildPdbFileStructureFromAssembly(-1,0);
	std::cout << "Writing output file:\n" << fullOutputFileName << "\n";
	outputPdbFile->Write(fullOutputFileName);
	delete outputPdbFile;
}

bool gmml::startsWith(std::string bigString, std::string smallString)
{
	return (bigString.compare(0, smallString.length(), smallString) == 0);
}

int gmml::CountInternalBonds(Assembly &ass)
{
	// This belong in a Molecule class. Algo:
	// Get heavy atoms in molecule. Check number of neighbors. Check if there are more atoms within distance d than neighbors. Scream if so.
	int count = 0;
	int withinDistance = 0;
	AtomVector heavyAtoms = selection::FindHeavyAtoms(ass.GetAllAtomsOfAssembly());
	for (auto &atom : heavyAtoms)
	{
		withinDistance = selection::AtomsWithinDistanceOf(atom, 1.9, heavyAtoms).size();
		//std::cerr << atom->GetId() << ": " << withinDistance << "\n";
		count += withinDistance;
	}
	return count;
}
