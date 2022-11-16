#include "includes/CentralDataStructure/atom.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp" // bondIfClose
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/MolecularMetadata/elementattributes.hpp"
#include <ctype.h> // isalpha

using cds::Atom;
//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Atom::Atom(const std::string& name, const Coordinate& coord)
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
    return this->getName();
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



 // This or just construct a pdbAtom class and write from there? No. This should be a free templated(?) function
void Atom::WritePdb(std::ostream &stream,
        std::string residueName,
        std::string residueNumber,
        std::string recordName,
        std::string chainId,
        std::string insertionCode,
        std::string alternativeLocation,
        std::string occupancy,
        std::string temperatureFactor
        ) const
{
    stream << std::left << std::setw(6) << recordName;
    stream << std::right << std::setw(5) << this->getIndex() << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(4) << this->getName();
    stream << std::left << std::setw(1) << alternativeLocation;
    stream << std::right << std::setw(3) << residueName << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(1) << chainId;
    stream << std::right << std::setw(4) << residueNumber;
    stream << std::left << std::setw(1) << insertionCode << std::left << std::setw(3) << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetX();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetY();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetZ();
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << occupancy;
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << temperatureFactor << std::left << std::setw(10) << " ";
    stream << std::right << std::setw(2) << this->getElement();
    if (this->getCharge() != codeUtils::dNotSet)
    {
        stream << std::left << std::setw(2) << std::to_string(this->getCharge()) << std::endl;
    }
    return;
}
