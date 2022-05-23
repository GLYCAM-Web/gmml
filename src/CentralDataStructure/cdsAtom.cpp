#include <ctype.h> // isalpha

#include "includes/CentralDataStructure/cdsAtom.hpp"
#include "includes/CentralDataStructure/cdsCoordinate.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp" // bondIfClose
#include "includes/CodeUtils/logging.hpp"

using cds::cdsAtom;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
cdsAtom::cdsAtom()
: Node<cdsAtom> ("DEFAULT CDS ATOM") {}

cdsAtom::cdsAtom(const std::string& name, const Coordinate& coord)
: Node<cdsAtom> (name)
{
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

const Coordinate* cdsAtom::getCoordinate() const
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
//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////

double cdsAtom::calculateDistance(const cdsAtom* otherAtom) const
{
    return this->getCoordinate()->Distance(otherAtom->getCoordinate());
}

void cdsAtom::addBond(cdsAtom* otherAtom)
{
    this->addNeighbor("bondByDistance", otherAtom);
}

std::string cdsAtom::getElement() // derived classes should overwrite if more explicit about element.
{
    std::string name = this->getName();
    if (!name.empty())
    {
        if (isalpha(name.at(0))) // if first char is in the alphabet
        {
             return name.substr(0,1); // return first character as string
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::WAR, "Did not find an element for atom named: " + name);
    return "";
}

void cdsAtom::bondIfClose(cdsAtom* otherAtom)
{
    double maxLength = atomicBonds::getMaxBondLengthByAtomType(this->getElement(), otherAtom->getElement());
    if (this->getCoordinate()->withinDistance(otherAtom->getCoordinate(), maxLength))
    {
        this->addBond(otherAtom);
        //std::stringstream ss;
        //std::cout << "Bonded " << this->getName() << "_" << this->getIndex() << " to " << otherAtom->getName() << "_" << otherAtom->getIndex() << "\n";
        //gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
    }
    return;
}
