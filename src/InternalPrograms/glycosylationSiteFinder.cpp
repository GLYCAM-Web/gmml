#include "includes/InternalPrograms/glycosylationSiteFinder.hpp"
#include "includes/MolecularModeling/residue.hpp"
#include "includes/MolecularMetadata/glycoprotein.hpp"

using glycoproteinBuilder::GlycosylationSiteFinder;

GlycosylationSiteFinder::GlycosylationSiteFinder(MolecularModeling::Assembly &ass)
{
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
        	table_.emplace_back(residue->GetChainID(), residue->GetNumber(), residue->GetInsertionCode(), context, tags );
    	}
    }
}

std::string GlycosylationSiteFinder::PrintTable()
{
    std::string output = "";
    for (auto &tableElement : table_)
    {
        output += tableElement.Print() + "\n";
    }
    return output;
}

