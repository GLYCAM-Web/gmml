#include "includes/MolecularMetadata/atomicBonds.hpp"
//#include "includes/CentralDataStructure/cdsAtom.hpp"

//#include <mutex>

using namespace atomicBonds;

bool atomicBonds::bondAtomsIfClose(cds::cdsAtom* atom1, cds::cdsAtom* atom2)
{
	//std::mutex mtx;           // mutex for critical section
    double maxLength = atomicBonds::getMaxBondLengthByAtomType(atom1->getElement(), atom2->getElement());
    if (atom1->getCoordinate()->withinDistance(atom2->getCoordinate(), maxLength))
    {
    	// std::lock_guard<std::mutex> guard(mtx); // pre C++17 version
    	//std::lock_guard guard(mtx);//RAII, the mutex will be unlocked upon guard destruction. Exception safe.
    	atom1->addBond(atom2);
    	//std::cout << "Bonded " << atom1->getName() << "_" << atom1->getIndex() << " to " << atom2->getName() << "_" << atom2->getIndex() << std::endl;
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
