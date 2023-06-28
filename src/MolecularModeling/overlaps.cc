#include "includes/MolecularModeling/overlaps.hpp"
#include "includes/MolecularModeling/atom.hpp"
#include "includes/MolecularModeling/atomnode.hpp"
#include "includes/common.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/constants.hpp"

using MolecularModeling::Atom;

double gmml::CalculateAtomicOverlaps(std::vector<Atom*> atomsA, std::vector<Atom*> atomsB, bool print)
{
    double totalOverlap   = 0.0;
    double currentOverlap = 0.0;
    for (auto& atomA : atomsA)
    {
        for (auto& atomB : atomsB)
        { // if not the same atom (index is unique)
            if ((atomA->GetIndex() != atomB->GetIndex()) && (atomA->CheckIfOtherAtomIsWithinOverlapDistance(atomB)))
            {
                currentOverlap = gmml::CalculateAtomicOverlaps(atomA, atomB);
                totalOverlap   += currentOverlap;
                if (print)
                {
                    //                	std::cout << atomA->GetId() << "::" << atomB->GetId() << ": " << (currentOverlap
                    //                / gmml::CARBON_SURFACE_AREA) << "\n";
                }
            }
        }
    }
    if (totalOverlap < 0.0) // Negative number fail
    {
        std::stringstream ss;
        ss << "Negative overlap should not happen, this is a bug: " << totalOverlap;
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return (totalOverlap / gmml::CARBON_SURFACE_AREA); // Normalise to area of a buried carbon
}

double gmml::CalculateAtomicOverlapsBetweenNonBondedAtoms(std::vector<Atom*>& atomsA, std::vector<Atom*>& atomsB)
{
    double totalOverlap = 0.0;
    for (auto& atomA : atomsA)
    {
        for (auto& atomB : atomsB)
        {
            bool isNeighbor              = false;
            std::vector<Atom*> neighbors = atomA->GetNode()->GetNodeNeighbors();
            for (auto& neighbor : neighbors)
            {
                if (atomB->GetIndex() == neighbor->GetIndex())
                {
                    isNeighbor = true;
                }
            }
            if ((isNeighbor == false) && (atomA->GetIndex() != atomB->GetIndex()) &&
                (atomA->CheckIfOtherAtomIsWithinBondingDistance(atomB)))
            {
                totalOverlap += gmml::CalculateAtomicOverlaps(atomA, atomB);
            }
        }
    }
    return (totalOverlap / gmml::CARBON_SURFACE_AREA); // Normalise to area of a buried carbon
}

double gmml::CalculateAtomicOverlaps(MolecularModeling::Atom* atomA, MolecularModeling::Atom* atomB, double radiusA,
                                     double radiusB)
{
    double distance = atomA->GetDistanceToAtom(atomB);
    if (radiusA == 0.0) // default value is 0.0, but user can provide.
    { // element info not usually set, so I look at first letter of atom name. This may be why you're reading this.
        if (atomA->GetName().at(0) == 'C')
        {
            radiusA = 1.70; // Rowland and Taylor modification to vdW.
        }
        if (atomA->GetName().at(0) == 'O')
        {
            radiusA = 1.52;
        }
        if (atomA->GetName().at(0) == 'N')
        {
            radiusA = 1.55;
        }
        if (atomA->GetName().at(0) == 'S')
        {
            radiusA = 1.80;
        }
        if (atomA->GetName().at(0) == 'P')
        {
            radiusA = 1.80;
        }
        if (atomA->GetName().at(0) == 'H')
        {
            radiusA = 1.09;
        }
    }
    if (radiusB == 0.0)
    {
        if (atomB->GetName().at(0) == 'C')
        {
            radiusB = 1.70;
        }
        if (atomB->GetName().at(0) == 'O')
        {
            radiusB = 1.52;
        }
        if (atomB->GetName().at(0) == 'N')
        {
            radiusB = 1.55;
        }
        if (atomB->GetName().at(0) == 'S')
        {
            radiusB = 1.80;
        }
        if (atomB->GetName().at(0) == 'P')
        {
            radiusB = 1.80;
        }
        if (atomB->GetName().at(0) == 'H')
        {
            radiusB = 1.09;
        }
    }
    // std::cout << "Distance: " << distance << " radiusA: " << radiusA << " radiusB: " << radiusB << std::endl;
    double overlap = 0.0;
    if (radiusA + radiusB > distance) // Close enough to overlap
    {
        if (std::abs(radiusA - radiusB) > distance) // If one sphere is completely inside the other
        {                                           // then calculate the surface area of the buried (smaller) sphere.
            if (radiusA < radiusB)
            {
                overlap = 4 * gmml::PI_RADIAN * (radiusA * radiusA);
            }
            else
            {
                overlap = 4 * gmml::PI_RADIAN * (radiusB * radiusB);
            }
        }
        else // Normal situation, partial overlap we need to calculate
        {    // Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours. Each atom against each atom, so
             // overlap can be "double" counted. See paper.
            overlap = (2 * (gmml::PI_RADIAN)*radiusA *
                       (radiusA - distance / 2 - (((radiusA * radiusA) - (radiusB * radiusB)) / (2 * distance))));
        }
    }
    if ((overlap < 0.0) || (radiusA == -0.1) || (radiusB == -0.1))
    { // Either the user didn't specify the radius or the element isn't one of the above
        //    	std::cout << "Neggie: " << overlap << " d: " << distance << ", A: " << atomA->GetName() << ", rA: " <<
        //    radiusA << ", B: " << atomB->GetName() << ", rB: " << radiusB << std::endl;
        return 0.0; // negative overlap isn't a thing.
    }
    // std::cout << "Non-normalized Overlap=" << totalOverlap << std::endl;
    return overlap;
}
