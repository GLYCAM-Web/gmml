#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_WIGGLETOSITE_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_WIGGLETOSITE_HPP_

#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/inputs.hpp"
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/assembly.hpp"

namespace gmmlPrograms
{
void wiggleToSite(WiggleToSiteInputs inputStruct);
void wiggleToSite(const std::vector<Coordinate*> atomsToAvoid, cdsCondensedSequence::Carbohydrate* carbohydrate, const cds::Residue* wiggleTarget, cds::Residue* wiggleMe, const cds::Residue* superimpositionTarget, cds::Residue* superimposeMe);
}
#endif /* INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_WIGGLETOSITE_HPP_ */
