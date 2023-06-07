#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_

#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/threads.hpp"
#include "includes/CodeUtils/constants.hpp" //dNotSet
#include "includes/CodeUtils/numbers.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <iomanip> // setprecision
#include <vector>
#include <thread>

namespace cds
{
void bondAtomsByDistance(std::vector<cds::Atom*> atoms);
double getCharge(std::vector<cds::Atom*> atoms);
void EnsureIntegralCharge(std::vector<cds::Atom*> atoms);
void bondAtomsByDistanceSerial(std::vector<cds::Atom*> atoms);
void bondAtomsAndResiduesByDistance(cds::Residue* residueA, cds::Residue* residueB);
void bondAtomsAndResiduesByDistance(std::vector<cds::Residue*> residues);
//Templated functions
template <typename T>
void serializeNumbers(std::vector<T*> elements)
{
    unsigned int i = 0;
    for(auto &element : elements)
    {
        element->setNumber(++i);
    }
    return;
}
} // namespace
#endif /* INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP_ */
