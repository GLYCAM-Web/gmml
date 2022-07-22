#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_

#include "includes/MolecularMetadata/atomicBonds.hpp"

#include <vector>
#include <thread> // thread


namespace cds
{
template <class atomT>
void bondAtomsByDistance(std::vector<atomT*> atoms)
{
//    std::vector<std::thread> threads;
//    int maxThreads = 11;
//    int i = 1;
	// number of threads from: std::thread::hardware_concurrency()
    for(typename std::vector<atomT*>::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        atomT* atom1 = *it1;
        for(typename std::vector<atomT*>::iterator it2 = std::next(it1); it2 != atoms.end(); ++it2)
        {
            atomT* atom2 = *it2;
            //gmml::log(__LINE__,__FILE__,gmml::INF, "Checking " + atom1->GetId() + " vs " + atom2->GetId());
            // Need to mutex lock the add node function in atom1? What about accessing atom2 as a node in that context? Probably also needs to lock.
            std::thread t(&atomicBonds::bondAtomsIfClose, atom1, atom2);
            //t.join(); This makes it linear.
            //atomicBonds::bondAtomsIfClose(atom1, atom2);
        }
    }
    return;
}
}



#endif /* INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_ */
