#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/coordinate.hpp"

using cds::Atom;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
Atom::Atom(const std::string& name, const Coordinate& coord)
{
    this->setName(name);
    this->addCoordinate(coord);
}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
cds::Coordinate* Atom::getCoordinate()
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
void Atom::setCoordinate(const Coordinate& newCoord)
{
    if(coordinates_.empty())
    {
        this->addCoordinate(newCoord);
    }
    else
    {
        Coordinate* firstCoord = coordinates_.front().get();
        firstCoord = newCoord; // Should copy the x,y,z from newCoord.
    }
    return;
}

void Atom::addCoordinate(const Coordinate& newCoord)
{
    coordinates_.push_back(std::make_unique<Coordinate>(newCoord));
    return;
}

