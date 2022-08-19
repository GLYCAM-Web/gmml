#ifndef INCLUDES_PARAMETERSET_CDS3DTEMPLATE_HPP
#define INCLUDES_PARAMETERSET_CDS3DTEMPLATE_HPP

#include "includes/CentralDataStructure/cdsResidue.hpp"
#include "includes/CentralDataStructure/cdsAtom.hpp"
#include "includes/CentralDataStructure/cdsFunctions.hpp"
#include "includes/ParameterSet/parameterManager.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"

// This may not end up being used in the new design. Leaving here just in case.

namespace cds
{

template <typename atomT>
std::vector<cdsResidue<atomT>> findAndCreateTemplatesForResidues(std::vector<std::string> queries)
{
	std::vector<cdsResidue<atomT>> results;
	parameters::Manager parameterManager; // loads in defaults.
	for (auto & query : queries)
	{
		PrepFileSpace::PrepFileResidue* prepResidue = parameterManager.FindPrepResidue(query);
		if (prepResidue == nullptr)
		{
			LibraryFileSpace::LibraryFileResidue* libResidue = parameterManager.FindLibResidue(query);
			if (libResidue == nullptr)
			{
				throw std::runtime_error("No lib or prep file entry for residue named: " + query);
			}
			results.push_back(cds::convertLibFileToResidue(libResidue));
		}
		else
		{
			results.push_back(convertPrepFileToResidue(prepResidue));
		}

	}
	return results;
}



}



#endif /* INCLUDES_PARAMETERSET_CDS3DTEMPLATE_HPP_ */
