#include "../../includes/gmml.hpp"
#include <iostream>
#include <string>
#include <algorithm> // find
#include <fstream>
#include <sstream>
#include <map>

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
};
//struct GlycosylationSiteTable
//{
//	std::vector<GlycosylationSiteInfo> oLinks_;
//	std::vector<GlycosylationSiteInfo> cLinks_;
//	std::vector<GlycosylationSiteInfo> nLinksAll_;
//	std::vector<GlycosylationSiteInfo> nLinksLikely_;
//	std::vector<GlycosylationSiteInfo> nLinksCysteine_;
//};
	// Functions
//bool isStringInVector(const std::string s, const std::vector<std::string> &v)
//{
//	if (std::find(v.begin(), v.end(), s) != v.end())
//		return true;
//	return false;
//}

std::string FindStringInStringMap(std::string s, std::unordered_map<std::string, std::string> sMap)
{
	std::unordered_map<std::string, std::string>::const_iterator found = sMap.find(s);
	if (found != sMap.end())
	{
		std::cout << "Found: " << found->first << " is a " << found->second << "\n";
		return found->second;
	}
	return "";
}

MolecularModeling::Residue* FindNeighborConnectedViaSpecificAtom(MolecularModeling::Residue *queryResidue, std::string queryAtomName)
{
	MolecularModeling::Residue* neighborResidue = nullptr;
	MolecularModeling::Atom* queryAtom = queryResidue->GetAtom(queryAtomName);
	if (queryAtom == nullptr)
	{
		std::cerr << "Error: An atom named " << queryAtomName << " not found in residue: " << queryResidue->GetId() << "\n";
		return nullptr;
	}
	for (auto &neighborAtom : queryAtom->GetNode()->GetNodeNeighbors())
	{
		if (neighborAtom->GetResidue() != queryResidue)
		{
			neighborResidue = neighborAtom->GetResidue();
		}
	}
	return neighborResidue;
}

std::string GetContextDetermineTags(MolecularModeling::Residue* residue, std::vector<std::string> &tags)
{
	std::string context = "";
	std::cout << "Neighbors (by residue connectivity) are:\n";
	MolecularModeling::Residue* firstPrecedingNeighbor = FindNeighborConnectedViaSpecificAtom(residue, "N");
	MolecularModeling::Residue* secondPrecedingNeighbor = nullptr;
	if (firstPrecedingNeighbor)
	{
		secondPrecedingNeighbor = FindNeighborConnectedViaSpecificAtom(firstPrecedingNeighbor, "N");
	}
	MolecularModeling::Residue* firstFollowingNeighbor = FindNeighborConnectedViaSpecificAtom(residue, "C");
	MolecularModeling::Residue* secondFollowingNeighbor = nullptr;
	if (firstFollowingNeighbor)
	{
		secondFollowingNeighbor = FindNeighborConnectedViaSpecificAtom(firstFollowingNeighbor, "C");
	}

	std::cout << "Context is:\n" << secondPrecedingNeighbor->GetId() << "\n" << firstPrecedingNeighbor->GetId() << "\n" << residue->GetId() << "\n"
			<< firstFollowingNeighbor->GetId() << "\n" << secondFollowingNeighbor->GetId() << "\n";
	return context;
}



//std::string findSequenceNeighborhood(std::vector<MolecularModeling::Residue*>::iterator it, std::vector<MolecularModeling::Residue*>::iterator begin, std::vector<MolecularModeling::Residue*>::iterator end)
//{
//	std::string context = "";
//	if (it != begin)
//	{ // safely go back one more
//		std::advance(it, -1);
//		if (it != begin)
//		{ // safely go back one more
//			std::advance(it, -1);
//			context += (*it)->GetName() + "_";
//			std::advance(it, 1);
//		}
//		context += (*it)->GetName() + "_";
//		std::advance(it, 1);
//	}
//	for(int i = 3; i > 0; i--)
//	{
//		if (it != end)
//		{
//			context += (*it)->GetName() + "_";
//			std::advance(it, 1);
//		}
//	}
//	return context;
//}

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

    std::unordered_map<std::string, std::string> aminoAcidNameToCodeMap ({
       												  	 { "ALA", "A" },
                                                         { "ARG", "R" },
                                                         { "ASN", "N" },
                                                         { "ASP", "D" },
                                                         { "CYS", "C" },
                                                         { "GLN", "Q" },
                                                         { "GLU", "E" },
                                                         { "GLY", "G" },
                                                         { "HIS", "H" },
                                                         { "ILE", "I" },
                                                         { "LEU", "L" },
                                                         { "LYS", "K" },
                                                         { "MET", "M" },
                                                         { "PHE", "F" },
                                                         { "PRO", "P" },
                                                         { "SER", "S" },
                                                         { "THR", "T" },
                                                         { "TRP", "W" },
                                                         { "TYR", "Y" },
                                                         { "VAL", "V" }, });

    std::unordered_map<std::string, std::string> residueLinkMap ({
    												  { "ASN", "nLink" },
                                                      { "THR", "oLink" },
                                                      { "SER", "oLink" },
                                                      { "TYR", "oLink" },
                                                      { "TRP", "cLink" },
                                                      { "NLN", "nLink" },
                                                      { "OLT", "oLink" },
                                                      { "OLS", "oLink" },
                                                      { "OLY", "oLink" },
                                                      { "CLW", "cLink" } });


    MolecularModeling::Assembly ass (argv[1], gmml::InputFileType::PDB);
    ass.BuildStructureByDistance();
    ass.GenerateResidueNodesInAssembly();
    std::vector<MolecularModeling::Residue*> residues = ass.GetResidues();
    for (std::vector<MolecularModeling::Residue*>::iterator it = residues.begin(); it != residues.end(); ++it)
    {

    	MolecularModeling::Residue* residue = (*it);
    	std::string linkType = FindStringInStringMap(residue->GetName(), residueLinkMap);
    	if (!linkType.empty())
    	{
        	std::cout << "Checking: " << residue->GetId() << "(" << residue->GetIndex() << ") with linkType: " << linkType << std::endl;

    		std::vector<std::string> tags;
        	std::string context = GetContextDetermineTags(residue, tags);
    		//tableInfo.emplace_back(residue->GetChain(), residue->GetInsertionCode(), residue->GetNumber(), sequenceContext, tags );
    	}
    }
    return 0;
}
