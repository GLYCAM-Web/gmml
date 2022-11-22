#include "includes/CentralDataStructure/Overlaps/cdsOverlaps.hpp"
#include "includes/CodeUtils/constants.hpp" // maxcutoff
using GeometryTopology::Coordinate;
// I meant to time which is faster
bool cds::CheckIfOtherCoordinateIsWithinDistance(const Coordinate* a, const Coordinate* b, const double distance)
{
    double xDiff = a->GetX() - b->GetX();
    double yDiff = a->GetY() - b->GetY();
    double zDiff = a->GetZ() - b->GetZ();
    if ( (xDiff*xDiff + yDiff*yDiff + zDiff*zDiff) < distance*distance )
    {
        return true;
    }
    return false;
}

bool cds::CheckIfOtherCoordinateIsWithinDistanceA(const Coordinate* a, const Coordinate* b, const double distance)
{
    double xDiff = std::abs(a->GetX() - b->GetX());
    double yDiff = std::abs(a->GetY() - b->GetY());
    double zDiff = std::abs(a->GetZ() - b->GetZ());
    if ( ( xDiff < distance ) && (yDiff < distance) && (zDiff < distance) )
    {
        return cds::CheckIfOtherCoordinateIsWithinDistance(a, b, distance);
    }
    return false;
}

double cds::CalculateAtomicOverlaps(cds::Atom* atomA, cds::Atom* atomB, double radiusA, double radiusB)
{
    double distance = atomA->getCoordinate()->Distance(atomB->getCoordinate());
    if (radiusA == 0.0) // default value is 0.0, but user can provide.
    {   // element info not usually set, so I look at first letter of atom name. This may be why you're reading this.
        if (atomA->getName().at(0) == 'C') radiusA = 1.70; // Rowland and Taylor modification to vdW.
        if (atomA->getName().at(0) == 'O') radiusA = 1.52;
        if (atomA->getName().at(0) == 'N') radiusA = 1.55;
        if (atomA->getName().at(0) == 'S') radiusA = 1.80;
        if (atomA->getName().at(0) == 'P') radiusA = 1.80;
        if (atomA->getName().at(0) == 'H') radiusA = 1.09;
    }
    if (radiusB == 0.0)
    {
        if (atomB->getName().at(0) == 'C') radiusB = 1.70;
        if (atomB->getName().at(0) == 'O') radiusB = 1.52;
        if (atomB->getName().at(0) == 'N') radiusB = 1.55;
        if (atomB->getName().at(0) == 'S') radiusB = 1.80;
        if (atomB->getName().at(0) == 'P') radiusB = 1.80;
        if (atomB->getName().at(0) == 'H') radiusB = 1.09;
    }
    //std::cout << "Distance: " << distance << " radiusA: " << radiusA << " radiusB: " << radiusB << std::endl;
    double overlap = 0.0;
    if (radiusA + radiusB > distance) // Close enough to overlap
    {
        if(std::abs(radiusA - radiusB) > distance) // If one sphere is completely inside the other
        { // then calculate the surface area of the buried (smaller) sphere.
            if (radiusA < radiusB)
            {
                overlap = 4 * constants::PI_RADIAN * (radiusA*radiusA);
            }
            else
            {
                overlap = 4 * constants::PI_RADIAN * (radiusB*radiusB);
            }
        }
        else // Normal situation, partial overlap we need to calculate
        { // Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours. Each atom against each atom, so overlap can be "double" counted. See paper.
            overlap = ( 2 * (constants::PI_RADIAN) * radiusA* ( radiusA - distance / 2 - ( ( (radiusA*radiusA) - (radiusB*radiusB) ) / (2 * distance) ) ) );
        }
    }
    if ( (overlap < 0.0) || (radiusA == -0.1) || (radiusB == -0.1) )
    { // Either the user didn't specify the radius or the element isn't one of the above
        std::cout << "Neggie: " << overlap << " d: " << distance << ", A: " << atomA->getName() << ", rA: " << radiusA << ", B: " << atomB->getName() << ", rB: " << radiusB << std::endl;
        return 0.0; // negative overlap isn't a thing.
    }
    //std::cout << "Non-normalized Overlap=" << totalOverlap << std::endl;
    return overlap;
}

double cds::CalculateAtomicOverlaps(std::vector<cds::Atom*> atomsA, std::vector<cds::Atom*> atomsB, bool print)
{
    double totalOverlap = 0.0;
    double currentOverlap = 0.0;
    for(auto &atomA : atomsA)
    {
        for(auto &atomB : atomsB)
        {   // if not the same atom (index is unique)
            if ( (atomA->getIndex() != atomB->getIndex()) && (cds::CheckIfOtherCoordinateIsWithinDistance(atomA->getCoordinate(), atomB->getCoordinate(), constants::maxCutOff*2)) )
            {
                currentOverlap = cds::CalculateAtomicOverlaps(atomA, atomB);
                totalOverlap += currentOverlap;
                if (print)
                {
                    std::cout << atomA->getId() << "::" << atomB->getId() << ": " << (currentOverlap / constants::CARBON_SURFACE_AREA) << "\n";
                }
            }
        }
    }
    if (totalOverlap < 0.0) // Negative number fail
    {
        std::stringstream ss;
        ss << "Negative overlap should not happen, this is a bug: " << totalOverlap;
        gmml::log(__LINE__,__FILE__,gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return (totalOverlap / constants::CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}
double cds::CalculateAtomicOverlapsBetweenNonBondedAtoms(std::vector<cds::Atom*>& atomsA, std::vector<cds::Atom*>& atomsB)
{
    double totalOverlap = 0.0;
    for(auto &atomA : atomsA)
    {
        for(auto &atomB : atomsB)
        {
            bool isNeighbor = false;
            std::vector<cds::Atom*> neighbors = atomA->getNeighbors();
            for(auto &neighbor : neighbors)
            {
                if (atomB->getIndex() == neighbor->getIndex())
                    isNeighbor = true;
            }
            if ( (isNeighbor == false) && (atomA->getIndex() != atomB->getIndex()) && (cds::CheckIfOtherCoordinateIsWithinDistance(atomA->getCoordinate(), atomB->getCoordinate(), constants::maxCutOff)))
            {
                totalOverlap += cds::CalculateAtomicOverlaps(atomA, atomB);
            }
        }
    }
    return (totalOverlap / constants::CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}

unsigned int cds::CountOverlappingAtoms(std::vector<cds::Atom*>& atomsA, std::vector<cds::Atom*>& atomsB)
{
    unsigned int overlappingAtoms = 1;
    for(auto &atomA : atomsA)
    {
        for(auto &atomB : atomsB)
        {
            if (cds::CheckIfOtherCoordinateIsWithinDistance(atomA->getCoordinate(), atomB->getCoordinate(), constants::maxCutOff))
            {
                ++overlappingAtoms;
            }
        }
    }
    return overlappingAtoms;
}
