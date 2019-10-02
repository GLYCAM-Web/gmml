#include "../../includes/MolecularModeling/overlaps.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/common.hpp"

double gmml::CalculateAtomicOverlaps(MolecularModeling::AtomVector atomsA, MolecularModeling::AtomVector atomsB){

    double distance = 0.0, totalOverlap = 0.0;
    for(auto &atomA : atomsA)
    {
        for(auto &atomB : atomsB)
        {
            if ( std::abs((atomA->GetCoordinates().at(0)->GetX() - atomB->GetCoordinates().at(0)->GetX())) < 3.6 ) // This is faster than calulating distance, and rules out tons of atom pairs.
            {
                distance = atomA->GetDistanceToAtom(atomB);
                if ( ( distance < 3.6 ) && ( distance > 0.0 ) ) //Close enough to overlap, but not the same atom
                {
//                    double overlap = gmml::CalculateAtomicOverlaps(atomA, atomB);
//                    if (overlap > 0.01)
//                        std::cout << atomA->GetId() << "--"  << atomB->GetId() << ": " << overlap << "\n";
//                    totalOverlap += overlap;
                    totalOverlap += gmml::CalculateAtomicOverlaps(atomA, atomB);
                }
            }
        }
    }
    return (totalOverlap / gmml::CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}


double gmml::CalculateAtomicOverlaps(MolecularModeling::Atom *atomA, MolecularModeling::Atom *atomB, double radiusA, double radiusB)
{
    double distance = atomA->GetDistanceToAtom(atomB);
    if (radiusA == -0.1) // default value is -0.1, but user can provide.
    {
        // element info not usually set, so I look at first letter of atom name. This may be why you're reading this.
        if (atomA->GetName().at(0) == 'C') radiusA = 1.70; // Rowland and Taylor modification to vdW.
        if (atomA->GetName().at(0) == 'O') radiusA = 1.52;
        if (atomA->GetName().at(0) == 'N') radiusA = 1.55;
        if (atomA->GetName().at(0) == 'S') radiusA = 1.80;
        if (atomA->GetName().at(0) == 'P') radiusA = 1.80;
        if (atomA->GetName().at(0) == 'H') radiusA = 1.09;
    }
    if (radiusB == -0.1) // default value is -0.1, but user can provide.
    {
        if (atomB->GetName().at(0) == 'C') radiusB = 1.70;
        if (atomB->GetName().at(0) == 'O') radiusB = 1.52;
        if (atomB->GetName().at(0) == 'N') radiusB = 1.55;
        if (atomB->GetName().at(0) == 'S') radiusB = 1.80;
        if (atomB->GetName().at(0) == 'P') radiusB = 1.80;
        if (atomB->GetName().at(0) == 'H') radiusB = 1.09;
    }
    //std::cout << "Distance: " << distance << " radiusA: " << radiusA << " radiusB: " << radiusB << std::endl;
    double overlap = 0.0;
    if (radiusA + radiusB > distance + 0.6){ // 0.6 overlap is deemed acceptable. (Copying chimera:)
        // Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours. Each atom against each atom, so overlap can be "double" counted. See paper.
        overlap = ( 2 * (gmml::PI_RADIAN) * radiusA* ( radiusA - distance / 2 - ( ( (radiusA*radiusA) - (radiusB*radiusB) ) / (2 * distance) ) ) );
    }
    //std::cout << "Non-normalized Overlap=" << totalOverlap << std::endl;
    return overlap;
}
