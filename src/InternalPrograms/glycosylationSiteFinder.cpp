#include "includes/InternalPrograms/glycosylationSiteFinder.hpp"
#include "includes/MolecularModeling/residue.hpp"
#include "includes/MolecularMetadata/glycoprotein.hpp"

using glycoproteinBuilder::GlycosylationSiteFinder;

GlycosylationSiteFinder::GlycosylationSiteFinder(std::vector<cds::Residue*> residues)
{
    for (auto &residue : residues)
    {
    	std::string linkType = glycoproteinMetadata::LookupLinkTypeForAminoAcidName(residue->getName());
    	if (!linkType.empty())
    	{
        	//std::cout << "Checking: " << residue->GetId() << "(" << residue->GetIndex() << ") with linkType: " << linkType << std::endl;
    		std::vector<std::string> tags = {linkType};
        	std::string context = glycoproteinMetadata::GetSequenceContextAndDetermineTags(residue, tags);
        	pdb::ResidueId residueId = residue->getId();
        	table_.emplace_back(residueId.getChainId(), residueId.getNumber(), residueId.getInsertionCode(), context, tags );
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

