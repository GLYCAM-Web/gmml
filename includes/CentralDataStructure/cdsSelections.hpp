#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_HPP

#include "includes/Abstract/absResidue.hpp"

#include <vector>

namespace cds
{
template< class residueT>
inline typename std::vector<residueT*> selectResiduesByType(std::vector<residueT*> inputResidues, Abstract::absResidue::Type queryType)
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

} // namespace

#endif // CDS SELECTIONS

