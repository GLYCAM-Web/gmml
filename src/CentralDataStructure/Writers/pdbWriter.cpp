#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include <iomanip> // setw

void cds::writeMoleculeToPdb(std::ostream& stream, const std::vector<cds::Residue*> residues)
{
    for (auto &residue : residues)
    {
        cds::writeResidueToPdb(stream, residue);
    }
    stream << "TER\n";
}

void cds::writeResidueToPdb(std::ostream& stream, const cds::Residue* residue, const std::string recordName, bool addTerCard)
{
    for(auto &atom : residue->getAtoms())
    {
        cds::writeAtomToPdb(stream, atom, recordName, residue->getName(), residue->getNumber());
    }
    if(addTerCard)
    {
        stream << "TER\n";
    }
}

// Used by the PdbAtom class and cds classes. Thus what it takes as input is a bit odd.
void cds::writeAtomToPdb(std::ostream& stream, const cds::Atom* atom, const std::string recordName, const std::string residueName, const int residueNumber, const std::string chainId, const std::string insertionCode, const double occupancy, const double temperatureFactor)
{
    std::string residueAlternativeLocation = ""; // If we ever need this to be anything else, change this function.
    stream << std::left << std::setw(6) << recordName;
    stream << std::right << std::setw(5) << atom->getIndex() << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(4) << atom->getName();
    stream << std::left << std::setw(1) << residueAlternativeLocation;
    stream << std::right << std::setw(3) << residueName << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(1) << chainId;
    stream << std::right << std::setw(4) << residueNumber;
    stream << std::left << std::setw(1) << insertionCode << std::left << std::setw(3) << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->getCoordinate()->GetX();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->getCoordinate()->GetY();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->getCoordinate()->GetZ();
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << occupancy;
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << temperatureFactor << std::left << std::setw(10) << " ";
    stream << std::right << std::setw(2) << atom->getElement();
//    We probably don't want to write charges into pdb file. Width allowed is 2.
//    if (atom->getCharge() != codeUtils::dNotSet)
//    {
//        stream << std::left << std::setw(2) << std::setprecision(1) << atom->getCharge();
//    }
    stream << std::endl;
    return;
}
