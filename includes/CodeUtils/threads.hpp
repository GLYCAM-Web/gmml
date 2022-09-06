#ifndef INCLUDES_CODEUTILS_THREADS_HPP_
#define INCLUDES_CODEUTILS_THREADS_HPP_

#include <thread>

// Notes for locking using mutexes:
// std::mutex mtx;
// std::lock_guard<std::mutex> guard(mtx); // pre C++17 version
// std::lock_guard guard(mtx);//RAII, the mutex will be unlocked upon guard destruction. Exception safe.

namespace codeUtils
{
inline unsigned long calculateNumberOfThreads(unsigned long const length)
{
	// Work out details of data and from that the number of threads to use
	//unsigned long const length = std::distance(first, last);
	if(!length)
	{ // If vector is empty
		return 0;
	}
	unsigned long const minTasksPerThread = 25;
	unsigned long const max_threads = (length + minTasksPerThread - 1) / minTasksPerThread;
	unsigned long hardware_threads = std::thread::hardware_concurrency();
	if (hardware_threads == 0)
	{
		hardware_threads = 4;
	}
	unsigned long num_threads = std::min(hardware_threads, max_threads); // Choose whatever is the lesser (of two weevils).
	return num_threads;
	//unsigned long const block_size = length / num_threads;
};
} // namespace
#endif /* INCLUDES_CODEUTILS_THREADS_HPP_ */
