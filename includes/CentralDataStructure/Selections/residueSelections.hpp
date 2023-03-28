#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSSELECTIONS_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include <vector>

//ToDo is this useful or would cdsResidue work fine?
namespace cdsSelections
{
using cds::Residue;
//template< class residueT>
//inline typename std::vector<Residue*> selectResiduesByType(std::vector<Residue*> inputResidues, cds::ResidueType queryType)
std::vector<Residue*> selectResiduesByType(std::vector<Residue*> inputResidues, cds::ResidueType queryType);
//{
//	std::vector<Residue*> selectedResidues;
//	for(auto & residue : inputResidues)
//	{
//		if ( residue->GetType() == queryType )
//		{
//			selectedResidues.push_back(residue);
//		}
//	}
//	return selectedResidues;
//}
//cds::Molecule* findMoleculeOfResidue(std::vector<cds::Molecule*> molecules, Residue* queryResidue);
unsigned int findHighestResidueNumber(std::vector<Residue*> residues);
} // namespace

#endif // CDS SELECTIONS

