#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include <iomanip> // setw

void cds::writeMoleculeToPdb(std::ostream& stream, const std::vector<cds::Residue*> residues)
{
    auto it = residues.begin();
    while(it != residues.end())
//    for (auto &residue : residues)
    {
        cds::writeResidueToPdb(stream, *it);
        if((*it)->GetType() != cds::ResidueType::Protein)
        {
            stream << "TER\n";
        }
        ++it;
    }
    --it;
    if((*it)->GetType() == cds::ResidueType::Protein)
    {
        stream << "TER\n";
    }
}

void cds::writeResidueToPdb(std::ostream& stream, const cds::Residue* residue, const std::string recordName)
{
    for(auto &atom : residue->getAtoms())
    {
        cds::writeAtomToPdb(stream, atom, recordName, residue->getName(), residue->getNumber());
    }
}

// Used by the PdbAtom class and cds classes. Thus what it takes as input is a bit odd.
void cds::writeAtomToPdb(std::ostream& stream, const cds::Atom* atom, const std::string recordName, const std::string residueName, const int residueNumber, const std::string chainId, const std::string insertionCode, const double occupancy, const double temperatureFactor)
{
    std::string residueAlternativeLocation = ""; // If we ever need this to be anything else, change this function.
    stream << std::left << std::setw(6) << recordName;
    stream << std::right << std::setw(5) << atom->getNumber() << std::left << std::setw(1) << " ";
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

void cds::writeConectCards(std::ostream& stream, std::vector<cds::Residue*> residues)
{ // These are only written for atoms connecting residues. The numbers overflow/truncate when longer than 5, but the format is what the format is.
    for(auto &residue : residues)
    {
        std::vector<std::pair<const Atom*,const Atom*>> atomsPairsConnectedToOtherResidues = residue->getAtomPairsConnectedToOtherResidues();
        for(auto &atomPair : atomsPairsConnectedToOtherResidues)
        {  // I hate that order matters, but here we are:
//            if (atomPair.first->getNumber() < atomPair.second->getNumber())
//            {
                stream << "CONECT" << std::right << std::setw(5) << atomPair.first->getNumber() << std::right << std::setw(5) << atomPair.second->getNumber() << "\n";
//            }
//            else
//            {
//                stream << "CONECT" << std::right << std::setw(5) << atomPair.second->getNumber() << std::right << std::setw(5) << atomPair.first->getNumber() << "\n";
//            }
        }
    }
}
