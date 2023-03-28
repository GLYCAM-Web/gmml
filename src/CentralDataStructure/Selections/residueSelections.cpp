#include "includes/CentralDataStructure/Selections/residueSelections.hpp"

using cds::Residue;

std::vector<Residue*> cdsSelections::selectResiduesByType(std::vector<Residue*> inputResidues, cds::ResidueType queryType)
{
    std::vector<Residue*> selectedResidues;
    for(auto & residue : inputResidues)
    {
        if ( residue->GetType() == queryType )
        {
            selectedResidues.push_back(residue);
        }
    }
    return selectedResidues;
}

//cds::Molecule* cdsSelections::findMoleculeOfResidue(std::vector<cds::Molecule*> molecules, Residue* queryResidue)
//{
//    for(auto & molecule : molecules)
//    {
//        std::vector<Residue*> residuesOfMolecule = molecule->getResidues();
//        if(std::find(residuesOfMolecule.begin(), residuesOfMolecule.end(), queryResidue) != residuesOfMolecule.end())
//        {
//            return molecule;
//        }
//    }
//    return nullptr;
//}

unsigned int cdsSelections::findHighestResidueNumber(std::vector<Residue*> residues)
{
    unsigned int highest = residues.back()->getNumber(); // Good start.
    for (auto &residue : residues)
    {
        unsigned int resNumInt = residue->getNumber();
        if (resNumInt > highest)
        {
            highest = resNumInt;
        }
    }
    return highest;
}
