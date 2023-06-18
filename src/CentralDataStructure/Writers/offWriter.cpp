#include "includes/CentralDataStructure/Writers/offWriter.hpp"

std::string cds::getOffType(const cds::ResidueType queryType)
{
    if (queryType == cds::ResidueType::Protein)
    {
        return "p";
    }
    if (queryType == cds::ResidueType::Solvent)
    {
        return "w";
    }
    return "?";
}

void cds::WriteOffFileUnit(std::vector<cds::Residue*> residues, std::ostream& stream, const std::string unitName)
{
    // WriteAtomSection
    const std::string FLAG = "131072";
    stream << "!entry." << unitName
           << ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg"
           << std::endl;
    for (auto& residue : residues)
    {
        unsigned int atomNumberInResidue = 1;
        for (auto& atom : residue->getAtoms())
        {
            stream << " \"" << atom->getName() << "\" "
                   << "\"" << atom->getType() << "\" "
                   << "0"
                   << " " << residue->getNumber() << " " << FLAG << " " << atomNumberInResidue << " "
                   << atom->getAtomicNumber() << " " << std::fixed << std::setprecision(6) << atom->getCharge()
                   << std::endl;
            atomNumberInResidue++;
        }
    }
    // WriteAtomPertInfoSection
    stream << "!entry." << unitName
           << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << std::endl;
    for (auto& residue : residues)
    {
        for (auto& atom : residue->getAtoms())
        {
            stream << " \"" << atom->getName() << "\" "
                   << "\"" << atom->getType() << "\" " << 0 << " " << -1 << " " << std::setprecision(1) << 0.0
                   << std::endl;
        }
    }
    // WriteBoundBoxSection
    stream << "!entry." << unitName << ".unit.boundbox array dbl" << std::endl;
    stream << " "
           << "-1.000000" << std::endl;
    stream << " " << 0.0 << std::endl;
    stream << " " << 0.0 << std::endl;
    stream << " " << 0.0 << std::endl;
    stream << " " << 0.0 << std::endl;
    // WriteChildSequenceSection
    stream << "!entry." << unitName << ".unit.childsequence single int" << std::endl;
    stream << " " << residues.size() + 1 << std::endl;
    // WriteConnectSection
    //  Note: this is silly but fine for most cases. If you're reading this it's because it mattered and you need to
    //  make it better.
    stream << "!entry." << unitName << ".unit.connect array int" << std::endl;
    stream << " " << 1 << std::endl;
    stream << " " << residues.back()->getAtoms().back()->getNumber() << std::endl;
    // WriteConnectivitySection
    stream << "!entry." << unitName << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;
    for (auto& residue : residues)
    {
        for (auto& atom : residue->getAtoms())
        {
            for (auto& neighbor : atom->getChildren())
            { // According to docs: (the *second* atom is the one with the larger index). So ordering
                if (atom->getNumber() > neighbor->getNumber())
                {
                    stream << " " << atom->getNumber() << " " << neighbor->getNumber() << " " << 1 << std::endl;
                }
                else
                {
                    stream << " " << neighbor->getNumber() << " " << atom->getNumber() << " " << 1 << std::endl;
                }
            }
        }
    }
    // WriteHierarchySection
    stream << "!entry." << unitName << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx"
           << std::endl;
    for (auto& residue : residues)
    {
        stream << " \""
               << "U"
               << "\""
               << " " << 0 << " "
               << "\""
               << "R"
               << "\""
               << " " << residue->getNumber() << std::endl;
        for (auto& atom : residue->getAtoms())
        {
            stream << " \""
                   << "R"
                   << "\""
                   << " " << residue->getNumber() << " "
                   << "\""
                   << "A"
                   << "\""
                   << " " << atom->getNumber() << std::endl;
        }
    }
    // WriteNameSection
    stream << "!entry." << unitName << ".unit.name single str" << std::endl;
    stream << " \"" << unitName << "\"" << std::endl;
    // WritePositionSection
    stream << "!entry." << unitName << ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;
    for (auto& residue : residues)
    {
        for (auto& atom : residue->getAtoms())
        {
            stream << std::setprecision(6) << std::fixed << " " << atom->getCoordinate()->GetX() << " "
                   << atom->getCoordinate()->GetY() << " " << atom->getCoordinate()->GetZ() << std::endl;
        }
    }
    // WriteResidueConnectSection // Every residue needs a head/tail regardless of reality. tleap uses this info.
    stream << "!entry." << unitName
           << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << std::endl;
    for (auto& residue : residues)
    {
        auto atomsConnectedToOtherResidues = residue->getAtomsConnectedToOtherResidues();
        // Deal with residues that don't have a tail/head in reality:
        if (atomsConnectedToOtherResidues.size() == 1)
        { // Repeating the same atom changes the tree structure in the parm7 file. Not sure anything uses that. Old gmml
          // code repeats so doing that.
            // For reducing terminal old code puts a 1 2 0 0 0 0. So not repeat. Changing to 2 2 0 0 0 0 causes first
            // atom not to be M (main). Might be an issue let's see.
            atomsConnectedToOtherResidues.push_back(atomsConnectedToOtherResidues.front());
        }
        for (auto& atom : atomsConnectedToOtherResidues)
        {
            stream << " " << atom->getNumber();
        }
        int columnsWithZero = 6 - atomsConnectedToOtherResidues.size();
        for (int i = 0; i < columnsWithZero; ++i)
        {
            stream << " "
                   << "0";
        }
        stream << std::endl;
    }
    // WriteResiduesSection
    stream << "!entry." << unitName
           << ".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx"
           << std::endl;
    for (auto& residue : residues)
    {
        unsigned int childseq   = residue->getAtoms().size() + 1;
        unsigned int startatomx = residue->getAtoms().front()->getNumber();
        std::string restype     = cds::getOffType(residue->GetType());
        unsigned int imagingx   = 0;
        stream << " \"" << residue->getName() << "\""
               << " " << residue->getNumber() << " " << childseq << " " << startatomx << " "
               << "\"" << restype << "\""
               << " " << imagingx << std::endl;
    }
    // WriteSolventCapSection
    stream << "!entry." << unitName << ".unit.solventcap array dbl" << std::endl;
    stream << " "
           << "-1.000000" << std::endl;
    stream << " "
           << "0.0" << std::endl;
    stream << " "
           << "0.0" << std::endl;
    stream << " "
           << "0.0" << std::endl;
    stream << " "
           << "0.0" << std::endl;
    // WriteVelocitiesSection
    stream << "!entry." << unitName << ".unit.velocities table  dbl x  dbl y  dbl z" << std::endl;
    for (auto& residue : residues)
    {
        std::vector<cds::Atom*> atoms = residue->getAtoms();
        for (std::vector<cds::Atom*>::iterator i = atoms.begin(); i != atoms.end(); ++i)
        { // Maybe later we'll want to deal with atom velocities...
            stream << " "
                   << "0.0"
                   << " "
                   << "0.0"
                   << " "
                   << "0.0" << std::endl;
        }
    }
    return;
}

// ToDo shouldn't this have a .lib suffix?
void cds::WriteResiduesToOffFile(std::vector<cds::Residue*> residues, std::ostream& stream)
{ // For writing each residue separately
    stream << "!!index array str" << std::endl;
    for (auto& residue : residues)
    {
        stream << " \"" << residue->getName() << "\"" << std::endl;
    }
    for (auto& residue : residues)
    {
        cds::serializeNumbers(std::vector<cds::Residue*> {residue});
        cds::serializeNumbers(residue->getAtoms());
        cds::WriteOffFileUnit(std::vector<cds::Residue*> {residue}, stream, residue->getName());
    }
    return;
}

void cds::WriteMoleculeToOffFile(cds::Molecule* molecule, std::ostream& stream, const std::string unitName)
{ // For writing residues together as a molecule
    stream << "!!index array str" << std::endl;
    stream << " \"" << unitName << "\"" << std::endl;
    cds::serializeNumbers(molecule->getAtoms());
    cds::serializeNumbers(molecule->getResidues());
    cds::WriteOffFileUnit(molecule->getResidues(), stream, unitName);
    return;
}

void cds::WriteAssemblyToOffFile(cds::Assembly* assembly, std::ostream& stream, const std::string unitName)
{
    stream << "!!index array str" << std::endl;
    stream << " \"" << unitName << "\"" << std::endl;
    cds::serializeNumbers(assembly->getAtoms());
    cds::serializeNumbers(assembly->getResidues());
    cds::WriteOffFileUnit(assembly->getResidues(), stream, unitName);
    return;
}
