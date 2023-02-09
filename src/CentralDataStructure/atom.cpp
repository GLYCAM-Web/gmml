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
: Abstract::absAtom(coord)
{
    this->setName(name);
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
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void Atom::addBond(Atom *otherAtom)
{
    this->addNeighbor("bondByDistance", otherAtom);
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
