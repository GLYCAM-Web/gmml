#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_

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

//ToDo drop the templates and switch to cds classes.
namespace cds
{
// Helper struct for next function.
// Haven't yet figured out why this needs to be a struct for std::thread. Probably it's to provide the type.
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
    gmml::log(__LINE__, __FILE__, gmml::INF, "Setting atom connectivity by distance.");
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
		if (atom->getCharge() != constants::dNotSet) // Better to throw/warn here?
		{
			totalCharge += atom->getCharge();
		}
	}
	return totalCharge;
}

template <typename atomT>
void EnsureIntegralCharge(std::vector<atomT*> atoms)
{
	double charge = cds::getCharge(atoms);
	std::stringstream ss;
	ss << std::fixed;
	ss << "Total charge is: " << std::setprecision(5) << charge << std::endl;
	gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
	if (!codeUtils::isNumberIntegral(charge))
	{
		std::stringstream errorMessage;
		errorMessage << "Non-integral charge (" << charge << "). You cannot run MD with this.\n";
		std::cerr << errorMessage.str();
		throw errorMessage.str();
	}
	return;
}

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



#endif /* INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_HPP_ */
