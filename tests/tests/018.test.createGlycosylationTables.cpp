#include <iostream>
#include <string>
#include <vector>
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/MolecularModeling/residue.hpp"
#include "includes/MolecularMetadata/glycoprotein.hpp"

// Structs
struct GlycosylationSiteInfo
{
	// Constructor
	GlycosylationSiteInfo(std::string chain, std::string residueNumber, std::string insertionCode, std::string sequenceContext, std::vector<std::string> tags) : chain_(chain), residueNumber_(residueNumber), insertionCode_(insertionCode), sequenceContext_(sequenceContext), tags_(tags) {}
	// Data
	std::string chain_;
	std::string residueNumber_;
	std::string insertionCode_;
	std::string sequenceContext_;
	std::vector<std::string> tags_; // e.g. oLink, nLink, sequon, cysteineSequon, all, etc
	// Print
	inline std::string Print()
	{
	    std::string output = "Chain: " + chain_ + "\nResidueNumber: " + residueNumber_ + "\nInsertionCode: " + insertionCode_ + "\nSequenceContext: " + sequenceContext_ + "\nTags:\n";
	    for (auto &tag : tags_ )
	    {
	        output += tag + "\n";
	    }
	    return output;
	}
};

// To become the test
int main(int argc, char* argv[])
{
    if ( (argc != 2) && (argc != 3) )
    {
        std::cout << "Usage: createGlycosylationTables.exe inputFile.pdb [outputFileName]\n";
        std::cout << "Example: pdb2glycam 1RVX.pdb \n";
        std::exit(EXIT_FAILURE);
    }
    std::vector<GlycosylationSiteInfo> tableInfo;
    MolecularModeling::Assembly ass (argv[1], gmml::InputFileType::PDB);
    ass.BuildStructureByDistance(10); // number of threads to use.
    ass.GenerateResidueNodesInAssembly();
    for (auto &residue : ass.GetResidues())
    {
    	std::string linkType = glycoproteinMetadata::LookupLinkTypeForAminoAcidName(residue->GetName());
    	if (!linkType.empty())
    	{
        	//std::cout << "Checking: " << residue->GetId() << "(" << residue->GetIndex() << ") with linkType: " << linkType << std::endl;
    		std::vector<std::string> tags = {linkType};
        	std::string context = glycoproteinMetadata::GetSequenceContextAndDetermineTags(residue, tags);
        	tableInfo.emplace_back(residue->GetChainID(), residue->GetNumber(), residue->GetInsertionCode(), context, tags );
    	}
    }
    for (auto &tableElement : tableInfo)
    {
        std::cout << tableElement.Print() << "\n";
    }
    return 0;
}
