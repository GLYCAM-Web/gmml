#ifndef INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP
#define INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP

#include <utility>
#include <map>

#include "includes/CentralDataStructure/cdsAtom.hpp"


namespace atomicBonds
{
const double maxCutOff = 1.65;
const double minCutOff = 0.7;

const std::map<std::string, std::pair<double, double> > bondLengthMap =
{
  {"CC", std::make_pair(1.22, 1.67)},
  {"CO", std::make_pair(1.07975, 1.67025)},
  {"OC", std::make_pair(1.07975, 1.67025)},
  {"CN", std::make_pair(1.26, 1.55)},
  {"NC", std::make_pair(1.26, 1.55)},
  {"OP", std::make_pair(1.35, 1.776)},
  {"PO", std::make_pair(1.35, 1.776)},
  {"OS", std::make_pair(1.43, 1.78)},
  {"SO", std::make_pair(1.43, 1.78)},
  {"NS", std::make_pair(1.62, 1.77)},
  {"SN", std::make_pair(1.62, 1.77)}
};
// FUNCTIONS
std::pair<double,double> getBondLengthByAtomType(const std::string& atom1Element, const std::string& atom2Element);
double getMaxBondLengthByAtomType(const std::string &atom1Element, const std::string &atom2Element);
bool bondAtomsIfClose(cds::cdsAtom* atom1, cds::cdsAtom* atom2);
}



#endif /* INCLUDES_MOLECULARMETADATA_ATOMICBONDS_HPP_ */
