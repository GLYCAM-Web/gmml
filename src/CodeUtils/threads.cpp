#include <thread>
#include "includes/CodeUtils/threads.hpp"

unsigned long codeUtils::calculateNumberOfThreads(const unsigned long length)
{
    // Work out details of data and from that the number of threads to use
    // unsigned long const length = std::distance(first, last);
    if (!length)
    { // If vector is empty
        return 0;
    }
    const unsigned long minTasksPerThread = 25;
    const unsigned long max_threads       = (length + minTasksPerThread - 1) / minTasksPerThread;
    unsigned long hardware_threads        = std::thread::hardware_concurrency();
    if (hardware_threads == 0)
    {
        hardware_threads = 4;
    }
    unsigned long num_threads =
        std::min(hardware_threads, max_threads); // Choose whatever is the lesser (of two weevils).
    return num_threads;
    // unsigned long const block_size = length / num_threads;
};
