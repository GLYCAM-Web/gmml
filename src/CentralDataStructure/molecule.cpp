#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
#include "includes/CentralDataStructure/Writers/cdsOffWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"

using cds::Molecule;
using cds::Residue;
using cds::Atom;
//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Molecule::Molecule(const std::string chainID) : number_(0), chainId_(chainID) {}
// Move Ctor
Molecule::Molecule(Molecule&& other) noexcept : glygraph::Node<cds::Molecule>(other)
{
    residues_ = std::move(other.residues_);
    number_ = std::move(other.number_);
}
// Copy Ctor
Molecule::Molecule(const Molecule& other) : glygraph::Node<cds::Molecule>(other),
number_(other.number_)
{
    for (auto& residue : other.residues_)
    {
        residues_.push_back(std::make_unique<Residue>((*residue.get())) );
    }
}
// Move and Copy assignment operator
Molecule& Molecule::operator=(Molecule other)
{
    glygraph::Node<cds::Molecule>::operator=(other); //ToDo ok?
    swap(*this, other);
    return *this;
}
//////////////////////////////////////////////////////////
//                    ACCESSOR                          //
//////////////////////////////////////////////////////////
std::vector<Residue*> Molecule::getResidues() const
{
    std::vector<Residue*> residues;
    for(auto &residuePtr : residues_)
    {
        residues.push_back(residuePtr.get());
    }
    return residues;
}

std::vector<Atom*> Molecule::getAtoms() const
{
    std::vector<Atom*> atoms;
    for(auto &residue : this->getResidues())
    {
        std::vector<Atom*> currentResidueAtoms = residue->getAtoms();
        // Concatenates the vectors. currentResidueAtoms isn't left in a defined state but that's ok here.
        atoms.insert( atoms.end(), std::make_move_iterator(currentResidueAtoms.begin()),
                std::make_move_iterator(currentResidueAtoms.end()) );
    }
    return atoms;
}
//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void Molecule::addResidue(std::unique_ptr<Residue> myResidue)
{ // This is good: myResidue contains a vector of unique_ptr, so you don't want to copy that.
    residues_.push_back(std::move(myResidue));
}

void Molecule::setResidues(std::vector<std::unique_ptr<Residue>> myResidues)
{ // This is good: myResidue contains a vector of unique_ptr, so you don't want to copy that.
    residues_ = std::move(myResidues);
}

Residue* Molecule::insertNewResidue(std::unique_ptr<Residue> myResidue, const Residue& positionReferenceResidue)
{
    auto position = this->findPositionOfResidue(&positionReferenceResidue);
    if (position != residues_.end())
    {
        gmml::log(__LINE__,__FILE__,gmml::INF, "New residue named " + myResidue->getName() + " will be born; You're welcome.");
        ++position; // it is ok to insert at end(). I checked. It was ok. Ok.
        position = residues_.insert(position, std::move(myResidue));
    }
    else
    {
        gmml::log(__LINE__,__FILE__,gmml::WAR, "Could not create residue named " + myResidue->getName() + " as referenceResidue was not found\n");
    }
    return (*position).get(); // Dereference the reference to a uniquePtr, then use get() to create a raw ptr...
}

std::vector<std::unique_ptr<Residue>>::iterator Molecule::findPositionOfResidue(const Residue* queryResidue)
{
    auto i = residues_.begin();
    auto e = residues_.end();
    while (i != e)
    {
        if (queryResidue == i->get())
        {
            return i;
        }
        else
        {
            ++i;
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find position of " + queryResidue->getName() + " in vector\n"); // every class should have a print?
    return e;
}

std::vector<Residue*> Molecule::getResidues(std::vector<std::string> queryNames)
{
    return codeUtils::getElementsWithNames(this->getResidues(), queryNames);
}

Residue* Molecule::getResidue(const std::string& queryName)
{
    return codeUtils::findElementWithName(this->getResidues(), queryName);
}

void Molecule::deleteResidue(Residue* residue)
{
    auto i = this->findPositionOfResidue(residue); // auto makes my life easier
    if (i != residues_.end())
    {
        gmml::log(__LINE__,__FILE__,gmml::INF, "Residue " + residue->getName() + " has been erased. You're welcome.");
        i = residues_.erase(i);
    }
    return;
}

void Molecule::renumberResidues(int newStartNumber)
{
    for(auto & residue : this->getResidues())
    {
        residue->setNumber(newStartNumber++);
    }
}

//////////////////////////////////////////////////////////
//                    DISPLAY                           //
//////////////////////////////////////////////////////////
void Molecule::WritePdb(std::ostream& stream) const
{
    cds::writeMoleculeToPdb(stream, this->getResidues());
}

void Molecule::WriteOff(std::ostream& stream) const
{
    cds::WriteMoleculeToOffFile(this, stream, this->getName());
}
