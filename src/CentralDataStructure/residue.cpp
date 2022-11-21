#include "includes/CentralDataStructure/residue.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/constants.hpp" // sNotSet
#include "includes/CodeUtils/templatedSelections.hpp"
#include "includes/CentralDataStructure/cdsFunctions.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
//#include "includes/CentralDataStructure/Shapers/geometryTopologyInterface.hpp"
#include "includes/InputSet/PdbFile/pdbResidueId.hpp"

using cds::Residue;
using cds::Atom;
//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
Residue::Residue(const std::string& residueName, const Residue *referenceResidue)
{
    this->setName(residueName);
    this->setNumber(referenceResidue->getNumber() + 1);
    this->determineType(residueName);
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

std::string Residue::getId(std::string moleculeNumber) const
{
    const std::string insertionCode = constants::sNotSet;
    pdb::ResidueId temp(this->getName(), std::to_string(this->getNumber()), insertionCode, moleculeNumber);
    return temp.print();
}

std::vector<Coordinate*> Residue::getCoordinates() const
{
    return cds::getCoordinatesFromAtoms(this->getAtoms());
}

Coordinate Residue::getGeometricCenter() const
{
    if (geometricCenter_.GetX() == codeUtils::dNotSet)
    {
        geometricCenter_ = this->calculateGeometricCenter();
    }
    return geometricCenter_;
}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////
void Residue::addAtom(std::unique_ptr<Atom> myAtom)
{
    atoms_.push_back(std::move(myAtom));
    return;
}

//void Residue::addAtom(Atom* myAtom)
//{
//    atoms_.push_back(std::make_unique<Atom>(myAtom));
//    return;
//}

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

bool Residue::contains(const Atom* queryAtom) const
{
    std::vector<Atom*> atoms = this->getAtoms();
    return ( std::find(atoms.begin(), atoms.end(), queryAtom) != atoms.end() );
}

std::vector<const Atom*> Residue::getAtomsConnectedToOtherResidues() const
{
    std::vector<const Atom*> foundAtoms;
    std::vector<Atom*> residueAtoms = this->getAtoms();
    for(auto &atom : residueAtoms)
    {
        for(auto &neighbor : atom->getNeighbors())
        { // check if neighbor is not one of the atoms in this residue.
            if(std::find(residueAtoms.begin(), residueAtoms.end(), neighbor) == residueAtoms.end())
            {
                foundAtoms.push_back(neighbor);
            }
        }
    }
    return foundAtoms;
}

void Residue::MakeDeoxy(std::string oxygenNumber)
{ // if oxygenNumber is 6, then C6-O6-H6O becomes C6-Hd
    Atom* hydrogenAtom = this->FindAtom("H" + oxygenNumber + "O");
    Atom* oxygenAtom = this->FindAtom("O" + oxygenNumber);
    Atom* carbonAtom = this->FindAtom("C" + oxygenNumber);
    // Add O and H charge to the C atom.
    carbonAtom->setCharge(carbonAtom->getCharge() + oxygenAtom->getCharge() + hydrogenAtom->getCharge());
    // Delete the H of O-H
    this->deleteAtom(hydrogenAtom);
    // Now transform the Oxygen to a Hd. Easier than deleting O and creating H. Note: this H looks weird in LiteMol as bond length is too long.
//    std::string newID = oxygenAtom->getId();
//    newID.replace(0,oxygenAtom->getName().size(),"Hd");
    oxygenAtom->setName("Hd");
    oxygenAtom->setType("H1");
    oxygenAtom->setCharge(0.0000);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Completed MakeDeoxy\n");
}

Coordinate Residue::calculateGeometricCenter() const
{
    return cds::calculateGeometricCenter(this->getCoordinates());
}
//////////////////////////////////////////////////////////
//                    DISPLAY                           //
//////////////////////////////////////////////////////////
void Residue::Print(std::ostream& out) const
{
    out << "Name: " << this->getName() << ", index:" << this->getIndex() << ":\n";
    for (auto &atom : this->getAtoms())
    {
        atom->Print(out);
        out << "\n";
    }
    std::cout << "\n";
}
