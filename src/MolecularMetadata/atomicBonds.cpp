#include "includes/MolecularMetadata/atomicBonds.hpp"
// #include "includes/CentralDataStructure/cdsAtom.hpp"

// #include <mutex>

using namespace atomicBonds;

std::pair<double, double> atomicBonds::getBondLengthByAtomType(const std::string& atom1Element,
                                                               const std::string& atom2Element)
{ // Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::string bothAtoms = atom1Element + atom2Element;

    if (bondLengthMap.find(bothAtoms) != bondLengthMap.end())
    {
        std::pair<double, double> cutoffDistances = bondLengthMap.at(bothAtoms);
        return cutoffDistances;
    }
    else
    {
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Using default binding cutoff of 1.65");
        return std::make_pair(atomicBonds::minCutOff, atomicBonds::maxCutOff);
    }
}

double atomicBonds::getMaxBondLengthByAtomType(const std::string& atom1Element, const std::string& atom2Element)
{ // Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::pair<double, double> result = atomicBonds::getBondLengthByAtomType(atom1Element, atom2Element);
    return result.second;
}
