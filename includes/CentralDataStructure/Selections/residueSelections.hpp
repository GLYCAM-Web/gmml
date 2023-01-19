#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSSELECTIONS_HPP

#include "includes/CentralDataStructure/residue.hpp"

#include <vector>

//ToDo drop the templates and switch to cds classes. ptr to cds::Residue or abs::Residue should be cool.
namespace cds
{
template< class residueT>
inline typename std::vector<residueT*> selectResiduesByType(std::vector<residueT*> inputResidues, cds::ResidueType queryType)
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

