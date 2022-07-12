#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_HPP

#include "includes/Abstract/absResidue.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>
namespace cds
{
template< class residueT>
inline typename std::vector<residueT*> selectProteinResidues(std::vector<residueT*> inputResidues, Abstract::absResidue::Type queryType)
{
	std::vector<residueT*> selectedResidues;
	for(auto & residue : inputResidues)
	{
		if ( residue->GetType() == queryType )
		{
			selectedResidues.push_back(residue);
		}
	}
	return selectedResidues;
}

//
//inline std::vector<Abstract::absResidue*> selectResidues(std::vector<Abstract::absResidue*> inputResidues, Abstract::absResidue::Type queryType) // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved")
//{
//	std::vector<Abstract::absResidue*> selectedResidues; // @suppress("Type cannot be resolved")
//	for(auto & residue : inputResidues)
//	{
//		if ( residue->GetType() == queryType )
//		{
//			selectedResidues.push_back(residue);
//		}
//	}
//	return selectedResidues;
//}
} // namespace

#endif // CDS SELECTIONS

