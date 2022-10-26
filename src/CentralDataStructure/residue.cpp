#include "includes/CentralDataStructure/residue.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
#include "includes/CentralDataStructure/cdsFunctions.hpp"

using cds::Residue;
using cds::Atom;
//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
Residue::Residue(const std::string& residueName, const Residue *referenceResidue)
{
    this->setName(residueName);
    this->setNumber(referenceResidue->getNumber() + 1);
}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////

std::vector<Atom*> Residue::getAtoms() const
{
    std::vector<Atom*> atoms;
    for(auto &atomPtr : atoms_)
    {
        atoms.push_back(atomPtr.get());
    }
    return atoms;
}

std::vector<std::string> Residue::getAtomNames() const
{
    std::vector<std::string> foundAtomNames;
    for(auto &atom : this->getAtoms())
    {
        foundAtomNames.push_back(atom->getName());
    }
    return foundAtomNames;
}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
void Residue::addAtom(std::unique_ptr<Atom> myAtom)
{
    atoms_.push_back(std::move(myAtom));
    return;
}

bool Residue::deleteAtom(const Atom* atom)
{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
    auto i = this->FindPositionOfAtom(atom); // auto makes my life easier
    if (i != atoms_.end())
    {
        gmml::log(__LINE__,__FILE__,gmml::INF, "Atom " + atom->getName() + " has been erased. You're welcome.");
        i = atoms_.erase(i); // this costs a lot as everything after i gets shifted.
        return true;
    }
    return false;
}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
typename std::vector<std::unique_ptr<Atom>>::iterator Residue::FindPositionOfAtom(const Atom* queryAtom)
{
    typename std::vector<std::unique_ptr<Atom>>::iterator i = atoms_.begin();
    typename std::vector<std::unique_ptr<Atom>>::iterator e = atoms_.end();
    while (i != e)
    {
        if (queryAtom == i->get())
        {
            return i;
        }
        else
        {
            ++i;
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find " + queryAtom->getName() + " in atom records\n");
    return e;
}

Atom* Residue::FindAtom(const std::string queryName) const
{
    return codeUtils::findElementWithName(this->getAtoms(), queryName);
}

Atom* Residue::FindAtom(const int& queryNumber) const
{
    return codeUtils::findElementWithNumber(this->getAtoms(), queryNumber);
}

std::vector<const Atom*> Residue::getAtomsConnectedToOtherResidues() const
{
    std::vector<const Atom*> foundAtoms;
    for(auto &atom : this->getAtoms())
    {
        for(auto &neighbor : atom->getNeighbors())
        { // check if neighbor is not one of the atoms in this residue.
            if(!codeUtils::isElementPresent(this->getAtoms().begin(), this->getAtoms().end(), neighbor) )
            {
                foundAtoms.push_back(neighbor);
            }
        }
    }
    return foundAtoms;
}
//////////////////////////////////////////////////////////
//                    DISPLAY                           //
//////////////////////////////////////////////////////////
void Residue::WritePdb(std::ostream& stream, const bool addTerCard) const
{
    for(auto &atom : this->getAtoms())
    {
        atom->WritePdb(stream);
    }
    if(addTerCard)
    {
        stream << "TER\n";
    }
}

void Residue::Print(std::ostream& out) const
{
    out << this->getName();
    for (auto &atom : this->getAtoms())
    {
        atom->Print(out);
    }
}
