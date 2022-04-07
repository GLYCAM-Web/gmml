#include "includes/CentralDataStructure/cdsAtom.hpp"
#include "includes/CentralDataStructure/cdsCoordinate.hpp"

using cds::cdsAtom;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
cdsAtom::cdsAtom(const std::string& name, const cdsCoordinate& coord)
{
    this->setName(name);
    this->addCoordinate(coord);
}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
cds::cdsCoordinate* cdsAtom::getCoordinate()
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
void cdsAtom::setCoordinate(const cdsCoordinate& newCoord)
{
    if(coordinates_.empty())
    {
        this->addCoordinate(newCoord);
    }
    else
    {
        cdsCoordinate* firstCoord = coordinates_.front().get();
        firstCoord = newCoord; // Should copy the x,y,z from newCoord.
    }
    return;
}

void cdsAtom::addCoordinate(const cdsCoordinate& newCoord)
{
    coordinates_.push_back(std::make_unique<cdsCoordinate>(newCoord));
    return;
}

