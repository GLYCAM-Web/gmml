#include "includes/CentralDataStructure/cdsAtom.hpp"
#include "includes/CentralDataStructure/cdsCoordinate.hpp"

using cds::cdsAtom;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
cdsAtom::cdsAtom(const std::string& name, const Coordinate& coord)
{
    this->setName(name);
    this->addCoordinate(coord);
}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
Coordinate* cdsAtom::getCoordinate()
{
    if(coordinates_.empty())
    {
        return nullptr;
    }
    return coordinates_.front().get();
}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
void cdsAtom::setCoordinate(const Coordinate& newCoord)
{ // Dealing with one coord, so want it to be the coord returned when someone calls getCoordiante.
    if(coordinates_.empty())
    {
        this->addCoordinate(newCoord);
    }
    else
    {
        Coordinate* firstCoord = coordinates_.front().get();
        *firstCoord = newCoord; // Should copy the x,y,z from newCoord.
    }
    return;
}

void cdsAtom::addCoordinate(const Coordinate& newCoord)
{
    coordinates_.push_back(std::make_unique<Coordinate>(newCoord));
    return;
}

