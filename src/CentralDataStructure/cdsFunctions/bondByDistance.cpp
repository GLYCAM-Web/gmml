#include <thread>
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp"
#include "includes/CodeUtils/threads.hpp"
#include "includes/CodeUtils/logging.hpp"

namespace cds
{ // Helper struct for next function.

    struct bondAtomsByDistanceThread
    {
        void operator()(std::vector<cds::Atom*>::iterator currentPosition, std::vector<cds::Atom*>::iterator currentEnd,
                        std::vector<cds::Atom*>::iterator compareStart, std::vector<cds::Atom*>::iterator compareEnd)
        { // Check every atom from current to end against every following atom.
          // std::vector<cds::Atom*>::iterator currentPosition = currentStart;
            unsigned long taskCounterTemp = 0;
            while (currentPosition != currentEnd)
            {
                cds::Atom* currentAtom                            = *currentPosition;
                std::vector<cds::Atom*>::iterator comparePosition = compareStart;
                while (comparePosition != compareEnd)
                {
                    cds::Atom* compareAtom = *comparePosition;
                    if (currentAtom != compareAtom) // if not the same atom.
                    {
                        static_cast<void>(
                            atomicBonds::bondAtomsIfClose(currentAtom, compareAtom)); // cast to ignore returned bool.
                    }
                    ++comparePosition;
                }
                ++currentPosition;
                ++taskCounterTemp;
            }
            //             std::cout << "Completed number of tasks: " << taskCounterTemp << std::endl;
            return;
        }
    };

    struct bondNthAtomByDistanceThread
    {
        void operator()(std::vector<cds::Atom*>::iterator currentPosition, const unsigned long jumpSize,
                        const unsigned long jumpsToDo, std::vector<cds::Atom*>::iterator compareEnd)
        { // Check every atom from current to end against every following atom.
          // std::vector<cds::Atom*>::iterator currentPosition = currentStart;
            unsigned long jumpCount      = 0;
            unsigned long atomsCompleted = 0;
            while (jumpCount < jumpsToDo)
            {
                cds::Atom* currentAtom                            = *currentPosition;
                std::vector<cds::Atom*>::iterator comparePosition = currentPosition + 1;
                while (comparePosition != compareEnd)
                {
                    cds::Atom* compareAtom = *comparePosition;
                    static_cast<void>(atomicBonds::bondAtomsIfClose(currentAtom, compareAtom));
                    ++comparePosition;
                }
                std::advance(currentPosition, jumpSize);
                ++atomsCompleted;
                ++jumpCount;
            }
            //            std::cout << "Atoms completed: " << atomsCompleted << std::endl;
            return;
        }
    };
} // namespace cds

void cds::bondAtomsByDistance(std::vector<cds::Atom*> atoms)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Setting atom connectivity by distance.");
    // Threading here by breaking up data into blocks.
    // Work out details of data and from that the number of threads to use
    const unsigned long vectorLength = atoms.size();
    // guesses from hardware, includes current load, not perfect.
    const unsigned long num_threads  = codeUtils::calculateNumberOfThreads(vectorLength);
    if (num_threads == 0)
    {
        const std::string message = "Tried to bondByDistance on an atom vector of length 0. Something is wrong.";
        gmml::log(__LINE__, __FILE__, gmml::INF, message);
        throw std::runtime_error(message);
        return;
    }
    // Break data up into blocks of data and launch threads
    std::vector<std::thread> threads(num_threads - 1);
    typename std::vector<cds::Atom*>::iterator block_start = atoms.begin();
    const unsigned long jumpsToDo                          = (vectorLength / num_threads);
    // std::cout << "vectorLength: " << vectorLength << ", num_threads: " << num_threads
    //            << ",jumpSize:" << num_threads - 1  << ", jumpstoDo: " << jumpsToDo << std::endl;

    for (unsigned long i = 0; i < (num_threads - 1); ++i)
    {
        //    std::cout << "Launching thread " << i << "\n";
        threads[i] =
            std::thread(cds::bondNthAtomByDistanceThread(), block_start, num_threads - 1, jumpsToDo, atoms.end());
        ++block_start;
    }
    // std::cout << "Completing the remaining tasks" << std::endl;
    //  Complete any remaining
    //     std::cout << "Will now advance " << jumpsToDo * (num_threads - 1) << std::endl;
    //     std::cout << "Remaining tasks should be " << vectorLength - (jumpsToDo * (num_threads - 1)) - (num_threads -
    //     1)
    //               << std::endl;
    std::advance(block_start, jumpsToDo * (num_threads - 1) - num_threads);
    cds::bondAtomsByDistanceThread()(block_start, atoms.end(), block_start, atoms.end());
    // Wait for all threads to finish before returning.
    for (auto& entry : threads)
    {
        entry.join();
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished setting atom connectivity by distance.");
    return;
}

// void cds::bondAtomsByDistance(std::vector<cds::Atom*> atoms)
//{
//     gmml::log(__LINE__, __FILE__, gmml::INF, "Setting atom connectivity by distance.");
//     // Threading here by breaking up data into blocks.
//     // Work out details of data and from that the number of threads to use
//     const unsigned long vectorLength = atoms.size();
//     // guesses from hardware, includes current load, not perfect.
//     const unsigned long num_threads = codeUtils::calculateNumberOfThreads(vectorLength);
//     if (num_threads == 0)
//     {
//         const std::string message = "Tried to bondByDistance on an atom vector of length 0. Something is wrong.";
//         gmml::log(__LINE__, __FILE__, gmml::INF, message);
//         throw std::runtime_error(message);
//         return;
//     }
//     const unsigned long block_size = vectorLength / num_threads;
//     // Break data up into blocks of data and launch threads
//     std::vector<std::thread> threads(num_threads - 1);
//     typename std::vector<cds::Atom*>::iterator block_start = atoms.begin();
//     for (unsigned long i = 0; i < (num_threads - 1); ++i)
//     {
//         // std::cout << "Launching thread " << i << "\n";
//         typename std::vector<cds::Atom*>::iterator block_end = block_start;
//         std::advance(block_end, block_size);
//         // std::cout << "Thread " << i << " checking from: " << (*block_start)->getIndex() << "_" <<
//         // (*block_start)->getNumber() << " to " << (*block_end)->getIndex() << std::endl;
//         threads[i]  = std::thread(cds::bondAtomsByDistanceThread(), block_start, block_end, block_start,
//         atoms.end()); block_start = block_end;
//     }
//     // Complete any remaining after division into blocks
//     cds::bondAtomsByDistanceThread()(block_start, atoms.end(), block_start, atoms.end());
//     // Wait for all threads to finish before returning.
//     for (auto& entry : threads)
//     {
//         entry.join();
//     }
//     gmml::log(__LINE__, __FILE__, gmml::INF, "Finished setting atom connectivity by distance.");
//     return;
// }

void cds::bondAtomsByDistanceSerial(std::vector<cds::Atom*> atoms)
{
    for (std::vector<cds::Atom*>::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        for (std::vector<cds::Atom*>::iterator it2 = std::next(it1); it2 != atoms.end(); ++it2)
        {
            atomicBonds::bondAtomsIfClose(*it1, *it2);
        }
    }
}

void cds::bondAtomsAndResiduesByDistance(cds::Residue* residueA, cds::Residue* residueB)
{
    bool residuesAreConnected = false;
    for (auto& atomA : residueA->getAtoms())
    {
        for (auto& atomB : residueB->getAtoms())
        {
            if (atomicBonds::bondAtomsIfClose(atomA, atomB))
            {
                residuesAreConnected = true; // only needs to be true once to connect residues.
            }
        }
    }
    if (residuesAreConnected)
    {
        std::string edgeName = residueA->getStringId() + "--" + residueB->getStringId();
        residueA->addNeighbor(edgeName, residueB);
    }
}

void cds::bondAtomsAndResiduesByDistance(std::vector<cds::Residue*> residues)
{
    for (std::vector<cds::Residue*>::iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
    { // First bond by distance for atoms within each residue
        cds::Residue* res1                = *it1;
        std::vector<cds::Atom*> res1Atoms = res1->getAtoms();
        cds::bondAtomsByDistanceSerial(res1Atoms);
        // Then for each residue, find other residues within reasonable residue distance.
        for (std::vector<cds::Residue*>::iterator it2 = std::next(it1); it2 != residues.end(); ++it2)
        {
            cds::Residue* res2                = *it2;
            std::vector<cds::Atom*> res2Atoms = res2->getAtoms();
            double residueDistance            = res1Atoms.at(0)->calculateDistance(res2Atoms.at(0));
            if (residueDistance < constants::residueDistanceOverlapCutoff * 1.2)
            {
                cds::bondAtomsAndResiduesByDistance(res1, res2);
            }
        }
    }
}

void cds::distanceBondIntra(std::vector<cds::Residue*> residues)
{
    for (auto& res : residues)
    {
        cds::bondAtomsByDistanceSerial(res->getAtoms());
    }
}

void cds::distanceBondInter(std::vector<cds::Residue*> residues)
{
    for (std::vector<cds::Residue*>::iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
    {
        cds::Residue* res1                = *it1;
        std::vector<cds::Atom*> res1Atoms = res1->getAtoms();
        for (std::vector<cds::Residue*>::iterator it2 = std::next(it1); it2 != residues.end(); ++it2)
        {
            cds::Residue* res2                = *it2;
            std::vector<cds::Atom*> res2Atoms = res2->getAtoms();
            double residueDistance            = res1Atoms.at(0)->calculateDistance(res2Atoms.at(0));
            if (residueDistance < constants::residueDistanceOverlapCutoff * 1.2)
            {
                cds::bondAtomsAndResiduesByDistance(res1, res2);
            }
        }
    }
}
