#include "includes/CentralDataStructure/atom.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp" // bondIfClose
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/MolecularMetadata/elementattributes.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include <ctype.h> // isalpha

using cds::Atom;
//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Atom::Atom(const std::string name, const Coordinate& coord)
{
    this->addCoordinate(coord);
    this->setName(name);
    this->setNumber(1); // Seems like a fine default?
}
// Move Ctor
Atom::Atom(Atom&& other) noexcept : glygraph::Node<cds::Atom>(other)
{
    coordinates_ = std::move(other.coordinates_);
    charge_ = std::move(other.charge_);
    atomType_ = std::move(other.atomType_);
    number_ = std::move(other.number_);
}
// Copy Ctor
Atom::Atom(const Atom& other) : glygraph::Node<cds::Atom>(other)
, charge_(other.charge_), atomType_(other.atomType_), number_(other.number_)
{
    for (auto& coord : other.coordinates_)
    {
        coordinates_.push_back( std::make_unique<Coordinate>((*coord.get())) );
    }
    //std::cout << "Atom ctor triggered\n";
}
// Move and Copy assignment operator
Atom& Atom::operator=(Atom other)
{
    glygraph::Node<cds::Atom>::operator=(other); //ToDo fuck.
    swap(*this, other);
    //std::cout << "MoveOrCopy operator triggered\n";
    return *this;
}

// This copy constructor causes const issues in Edge<T>
//Atom::Atom(Atom* refAtom) : Abstract::absAtom(refAtom->getCoordinate()), Node(*refAtom)
//{
//    this->setName(refAtom->getName());
//    this->setType(refAtom->getType());
//    this->setCharge(refAtom->getCharge());
//    this->setNumber(refAtom->getNumber());
//}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
// ToDo why is there two? const ref is better?
Coordinate* Atom::getCoordinate()
{
    if(coordinates_.empty())
    {
        return nullptr;
    }
    return coordinates_.front().get();
}

Coordinate* Atom::getCoordinate() const
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

void Atom::addCoordinate(const Coordinate& newCoord)
{
    coordinates_.push_back(std::make_unique<Coordinate>(newCoord));
    return;
}
//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void Atom::addBond(Atom *otherAtom)
{
    this->addNeighbor("atomicBond", otherAtom);
    return;
}

void Atom::bondIfClose(Atom* otherAtom)
{
    if (this->isWithinBondingDistance(otherAtom))
    {
        this->addBond(otherAtom);
        //std::stringstream ss;
        //std::cout << "Bonded " << this->getName() << "_" << this->getIndex() << " to " << otherAtom->getName() << "_" << otherAtom->getIndex() << "\n";
        //gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
    }
    return;
}

std::string Atom::getElement() const // derived classes should overwrite if more explicit about element.
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

int Atom::getAtomicNumber() const
{
    return MolecularMetadata::findElementAtomicNumber(this->getElement());
}

std::string Atom::getId() const
{
    return this->getName() + "_" + std::to_string(this->getIndex());
}

bool Atom::isWithinBondingDistance(const Atom* otherAtom) const
{
    double maxLength = atomicBonds::getMaxBondLengthByAtomType(this->getElement(), otherAtom->getElement());
    if (this->getCoordinate()->withinDistance(otherAtom->getCoordinate(), maxLength))
    {
        return true;
    }
    return false;
}

double Atom::calculateDistance(const Atom* otherAtom) const
{
    return this->getCoordinate()->Distance(otherAtom->getCoordinate());
}
//////////////////////////////////////////////////////////
//                   OVERLOADED OPERATORS               //
//////////////////////////////////////////////////////////
bool Atom::operator== (const Atom &otherAtom)
{
    return (this->getIndex() == otherAtom.getIndex());
}

bool Atom::operator!= (const Atom &otherAtom)
{
    return (this->getIndex() != otherAtom.getIndex());
}
//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void Atom::Print(std::ostream& out) const
{
    out << this->getName() << ", ";
    return;
}
