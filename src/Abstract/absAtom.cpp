#include "includes/Abstract/absAtom.hpp"

using Abstract::absAtom;
//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
absAtom::absAtom(const Coordinate& coord)
{
    this->addCoordinate(coord);
}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
// ToDo why is there two? const ref is better?
Coordinate* absAtom::getCoordinate()
{
    if(coordinates_.empty())
    {
        return nullptr;
    }
    return coordinates_.front().get();
}

const Coordinate* absAtom::getCoordinate() const
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
void absAtom::setCoordinate(const Coordinate& newCoord)
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

void absAtom::addCoordinate(const Coordinate& newCoord)
{
    coordinates_.push_back(std::make_unique<Coordinate>(newCoord));
    return;
}
//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
double absAtom::calculateDistance(const absAtom* otherAtom) const
{
    return this->getCoordinate()->Distance(*(otherAtom->getCoordinate()));
}


