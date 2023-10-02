#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp"

using pdb::residueSelector;

pdb::PdbResidue* pdb::residueSelector(const PdbFile& pdbFile, const pdb::ResidueId& residueId, const int modelNumber)
{
    for (auto& model : pdbFile.getAssemblies())
    {
        // std::cout << model->getNumber() << ":" << modelNumber << std::endl;
        if (model->getNumber() == modelNumber)
        {
            return pdb::residueSelector(model->getResidues(), residueId);
        }
    }
    return nullptr;
}

pdb::PdbResidue* pdb::residueSelector(std::vector<cds::Residue*> residues, const pdb::ResidueId& queryId)
{ // I'm using empty() to mean that it could be anything.
    for (auto& residue : residues)
    {
        PdbResidue* pdbResidue = static_cast<PdbResidue*>(residue);
        // std::cout << "currentId vs queryId: " << pdbResidue->getId() << " vs " << queryId << std::endl;
        if (queryId.getName().empty() || queryId.getName() == pdbResidue->getId().getName())
        {
            if (queryId.getNumber().empty() || queryId.getNumber() == pdbResidue->getId().getNumber())
            {
                if (queryId.getInsertionCode().empty() ||
                    queryId.getInsertionCode() == pdbResidue->getId().getInsertionCode())
                {
                    if ((queryId.getChainId().empty() || queryId.getChainId() == pdbResidue->getId().getChainId()))
                    {
                        return pdbResidue;
                    }
                }
            }
        }
    }
    return nullptr;
}
