#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_

#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/threads.hpp"

#include <vector>
#include <thread>

namespace cds
{
// Haven't yet figured out why this needs to be a struct for std::thread.
template <typename RandomAccessIterator, typename atomT> // atom1 passed in to give a type to the compiler
struct bondAtomsByDistanceThread
{
void operator() (RandomAccessIterator current, RandomAccessIterator end)
{ // Check every atom from current to end against every following atom.
	while (current != end)
	{
		atomT* atom1 = *current;
		for(typename std::vector<atomT*>::iterator it2 = std::next(current); it2 != end; ++it2)
		{
			atomT* atom2 = *it2;
			atomicBonds::bondAtomsIfClose(atom1, atom2);
		}
		++current;
	}
	return;
}
};

// ToDo exception handling
// Goal is to break data up into blocks and have each thread handle a block
template <typename atomT>
void bondAtomsByDistance(std::vector<atomT*> atoms)
{
	// Threading here by breaking up data into blocks.
	typename std::vector<atomT*>::iterator first = atoms.begin();
	typename std::vector<atomT*>::iterator last = atoms.end();
	// Work out details of data and from that the number of threads to use
	unsigned long const length = std::distance(first, last);
	unsigned long const num_threads = codeUtils::calculateNumberOfThreads(length);
	if (num_threads == 0) // if length is 0 (or 1?)
	{
		return;
	}
	unsigned long const block_size = length / num_threads;
	// Break data up into blocks of data and launch threads
	std::vector<std::thread> threads(num_threads - 1);
	typename std::vector<atomT*>::iterator block_start = first;
	for(unsigned long i = 0; i < (num_threads-1); ++i)
	{
		std::cout << "Launching thread " << i << "\n";
		typename std::vector<atomT*>::iterator block_end = block_start;
		std::advance(block_end, block_size);
		threads[i] = std::thread(cds::bondAtomsByDistanceThread<typename std::vector<atomT*>::iterator, atomT>(), block_start, block_end);
		block_start = block_end;
	}
	// Complete any remaining after division into blocks
	cds::bondAtomsByDistanceThread<typename std::vector<atomT*>::iterator, atomT>()(block_start, last);
	// Wait for all threads to finish before returning.
	std::cout << "Waiting on threads to finish\n";
    for(auto& entry: threads)
    {
    	entry.join();
    }
	std::cout << "Threads finish\n";
    return;
}

template <typename atomT>
double getCharge(std::vector<atomT*> atoms)
{
	double totalCharge = 0.0;
	for(auto &atom : atoms)
	{
		if (atom->getCharge() != gmml::dNotSet)
		{
			totalCharge += atom->getCharge();
		}
	}
	return totalCharge;
}

} // namespace



#endif /* INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_ */
