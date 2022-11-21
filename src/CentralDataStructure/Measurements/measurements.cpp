#include "includes/CentralDataStructure/Measurements/measurements.hpp"

Coordinate cds::calculateGeometricCenter(const std::vector<Coordinate*> coords)
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    for(auto &coord : coords)
    {
        x += coord->GetX();
        y += coord->GetY();
        z += coord->GetZ();
    }
    x = x / coords.size();
    y = y / coords.size();
    z = z / coords.size();
    return Coordinate(x,y,z);
}

std::vector<Coordinate*> cds::getCoordinatesFromAtoms(std::vector<cds::Atom*> atoms)
{
    std::vector<Coordinate*> coordinates;
    for(auto & atom : atoms)
    {
        coordinates.push_back(atom->getCoordinate());
    }
    return coordinates;
}
