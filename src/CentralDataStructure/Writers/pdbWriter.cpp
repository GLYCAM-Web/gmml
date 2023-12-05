#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include <iomanip> // setw

void cds::writeEnsembleToPdb(std::ostream& stream, const std::vector<cds::Assembly*> assemblies)
{
    int modelCount = 1;
    for (auto& assembly : assemblies)
    {
        stream << "MODEL " << std::right << std::setw(8) << modelCount << "\n";
        cds::writeAssemblyToPdb(stream, assembly->getMolecules());
        stream << "ENDMDL\n";
        modelCount++;
    }
}

void cds::writeAssemblyToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules)
{
    for (auto& molecule : molecules)
    {
        cds::writeMoleculeToPdb(stream, molecule->getResidues());
    }
}

void cds::writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules)
{
    unsigned int coordinateSets = molecules.at(0)->getAtoms().at(0)->getNumberOfCoordinateSets();
    for (unsigned int modelCount = 1; modelCount <= coordinateSets; modelCount++)
    {
        stream << "MODEL " << std::right << std::setw(8) << modelCount << "\n";
        for (auto& molecule : molecules)
        {
            cds::writeMoleculeToPdb(stream, molecule->getResidues(), (modelCount - 1));
        }
        stream << "ENDMDL\n";
    }
}

void cds::writeMoleculeToPdb(std::ostream& stream, const std::vector<cds::Residue*> residues,
                             unsigned int coordinateSetNumber)
{
    auto it = residues.begin();
    while (it != residues.end())
    {
        cds::writeResidueToPdb(stream, *it, "ATOM", coordinateSetNumber);
        if ((++it != residues.end()) && ((*it)->GetType() != cds::ResidueType::Protein))
        {
            stream << "TER\n";
        }
    }
    stream << "TER\n";
    //    --it;
    //    if((*it)->GetType() == cds::ResidueType::Protein)
    //    {
    //        stream << "TER\n";
    //    }
}

void cds::writeResidueToPdb(std::ostream& stream, const cds::Residue* residue, const std::string recordName,
                            unsigned int coordinateSetNumber)
{
    for (auto& atom : residue->getAtoms())
    {
        atom->setCurrentCoordinate(coordinateSetNumber);
        cds::writeAtomToPdb(stream, atom, recordName, residue->getName(), residue->getNumber());
    }
}

// Used by the PdbAtom class and cds classes. Thus what it takes as input is a bit odd.
void cds::writeAtomToPdb(std::ostream& stream, const cds::Atom* atom, const std::string recordName,
                         const std::string residueName, const int residueNumber, const std::string chainId,
                         const std::string insertionCode, const double occupancy, const double temperatureFactor)
{
    std::string residueAlternativeLocation = ""; // If we ever need this to be anything else, change this function.
    stream << std::left << std::setw(6) << recordName;
    stream << std::right << std::setw(5) << atom->getNumber() << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(4) << atom->getName();
    stream << std::left << std::setw(1) << residueAlternativeLocation;
    stream << std::right << std::setw(3) << residueName.substr(0, 3) << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(1) << chainId;
    stream << std::right << std::setw(4) << residueNumber;
    stream << std::left << std::setw(1) << insertionCode << std::left << std::setw(3) << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->getCoordinate()->GetX();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->getCoordinate()->GetY();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->getCoordinate()->GetZ();
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << occupancy;
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << temperatureFactor << std::left
           << std::setw(10) << " ";
    stream << std::right << std::setw(2) << atom->getElement();
    //    We probably don't want to write charges into pdb file. Width allowed is 2.
    //    if (atom->getCharge() != codeUtils::dNotSet)
    //    {
    //        stream << std::left << std::setw(2) << std::setprecision(1) << atom->getCharge();
    //    }
    stream << std::endl;
    return;
}

void cds::writeConectCards(std::ostream& stream, const std::vector<cds::Residue*> residues)
{ // These are only written for atoms connecting residues. The numbers overflow/truncate when longer than 5, but the
  // format is what the format is.
    std::vector<std::pair<const Atom*, const Atom*>> atomsPairsConnectedToOtherResidues;
    for (auto& residue : residues)
    {
        residue->findAtomPairsConnectedToOtherResidues(atomsPairsConnectedToOtherResidues);
    }
    //    auto it = std::unique(atomsPairsConnectedToOtherResidues.begin(),atomsPairsConnectedToOtherResidues.end() );
    //    atomsPairsConnectedToOtherResidues.resize(std::distance(atomsPairsConnectedToOtherResidues.begin(), it));
    for (auto& atomPair : atomsPairsConnectedToOtherResidues)
    {
        stream << "CONECT" << std::right << std::setw(5) << atomPair.first->getNumber() << std::right << std::setw(5)
               << atomPair.second->getNumber() << "\n";
        stream << "CONECT" << std::right << std::setw(5) << atomPair.second->getNumber() << std::right << std::setw(5)
               << atomPair.first->getNumber() << "\n";
    }
}
