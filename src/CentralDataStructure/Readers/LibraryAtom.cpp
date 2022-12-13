#include "includes/CentralDataStructure/Readers/LibraryAtom.hpp"
#include "includes/CodeUtils/strings.hpp"
using lib::LibraryAtom;

LibraryAtom::LibraryAtom (const std::string& line)
{
    int int_t; // Do not save this value in the class.
    std::string name, type;
    double charge;
    std::stringstream ss(line);
    ss >> name >> type >> int_t >> residue_index_ >> int_t >> atom_index_ >> atomic_number_ >> charge;     /// Split the line by space to extract attributes of the atom
    codeUtils::RemoveQuotes(name);
    codeUtils::RemoveQuotes(type);
    this->setName(name);
    this->setType(type);
    this->setCharge(charge);
    this->setNumber(this->getAtomIndex()); // Lib files' index is gmml's number.
    return;
}

int LibraryAtom::getResidueIndex()
{
    return residue_index_;
}

int LibraryAtom::getAtomIndex()
{
    return atom_index_;
}

int LibraryAtom::getAtomicNumber()
{
    return atomic_number_;
}

