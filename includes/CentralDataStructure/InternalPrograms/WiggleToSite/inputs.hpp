#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_INPUTS_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_INPUTS_HPP

#include <string>

namespace gmmlPrograms
{
struct WiggleToSiteInputs
{
    //ctor
    WiggleToSiteInputs(std::string inputFileName);
    //Members
	std::string carbohydrateSequence_ = "";
	int carbohydrateSuperimpositionResidue_ = 0;
	int carbohydrateWigglingResidue_ = 0;
	std::string substrateFile_ = "";
	int superimpositionTargetResidue_ = 0;
	int wigglingTargetResidue_ = 0;
	int persistCycles_ = 5;
	bool isDeterministic_ = false;
	std::string Print();
};
} // namespace
#endif // GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_WIGGLETOSITE_INPUTS_HPP
