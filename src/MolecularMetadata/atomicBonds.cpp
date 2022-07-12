#include "includes/MolecularMetadata/atomicBonds.hpp"
//#include "includes/CentralDataStructure/cdsAtom.hpp"

using namespace atomicBonds;

bool atomicBonds::bondAtomsIfClose(cds::cdsAtom* atom1, cds::cdsAtom* atom2)
{
    double maxLength = atomicBonds::getMaxBondLengthByAtomType(atom1->getElement(), atom2->getElement());
    if (atom1->getCoordinate()->withinDistance(atom2->getCoordinate(), maxLength))
    {
        atom1->addBond(atom2);
        //std::stringstream ss;
        //std::cout << "Bonded " << atom1->getName() << "_" << atom1->getIndex() << " to " << atom2->getName() << "_" << atom2->getIndex() << "\n";
        //gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
    }
    return false;
}

std::pair<double,double> atomicBonds::getBondLengthByAtomType(const std::string& atom1Element, const std::string& atom2Element)
{//Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::string bothAtoms = atom1Element + atom2Element;

    if(bondLengthMap.find(bothAtoms) != bondLengthMap.end())
    {
        std::pair<double,double> cutoffDistances = bondLengthMap.at(bothAtoms);
        return cutoffDistances;
    }
    else
    {
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Using default binding cutoff of 1.65");
        return std::make_pair(atomicBonds::minCutOff, atomicBonds::maxCutOff);
    }
}

double atomicBonds::getMaxBondLengthByAtomType(const std::string &atom1Element, const std::string &atom2Element)
{//Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::pair<double,double> result = atomicBonds::getBondLengthByAtomType(atom1Element, atom2Element);
    return result.second;
}
