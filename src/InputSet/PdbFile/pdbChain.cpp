#include "includes/InputSet/PdbFile/pdbChain.hpp"
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp" // get_cartesian_point_from_internal_coords
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/biology.hpp" // proteinResidueNames

using pdb::PdbChain;
////////////////////////////////////////////////////////////
////                       CONSTRUCTOR                    //
////////////////////////////////////////////////////////////
PdbChain::PdbChain(std::stringstream &stream_block, const std::string& chainId) : chainId_(chainId)
{
    gmml::log(__LINE__,__FILE__,gmml::INF, "Constructing PdbChain from stream_block >>>>>>>>>:\n" + stream_block.str() + "\n<<<<<<<<<<<<<< end stream_block");
    std::string line;
    while(getline(stream_block, line))
    {
        //gmml::log(__LINE__,__FILE__,gmml::INF, "Constructing chain with line: " + line);
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        if ( (recordName == "ATOM") || (recordName == "HETATM") )
        {
            std::stringstream singleResidueSection = this->extractSingleResidueFromRecordSection(stream_block, line);
            this->addResidue(std::make_unique<PdbResidue>(singleResidueSection, line));
        }
        else
        {
            gmml::log(__LINE__,__FILE__,gmml::INF, "In PdbChain Constructor with record that isn't cool: " + recordName);
            break;
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::INF, "Adding NTerminal and CTerminal tags if protein present");
    this->tagTerminalResidues();
    gmml::log(__LINE__,__FILE__,gmml::INF, "PdbChain Constructor Complete Captain");
    return;
}

////////////////////////////////////////////////////////////
////                         ACCESSOR                     //
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////                    FUNCTIONS                         //
////////////////////////////////////////////////////////////

void PdbChain::tagTerminalResidues()
{
    PdbResidue* nTer = this->getNTerminal();
    if (nTer != nullptr)
    {
        nTer->addLabel("NTerminal");
    }
    PdbResidue* cTer = this->getCTerminal();
    if (cTer != nullptr)
    {
        cTer->addLabel("CTerminal");
    }
    return;
}

std::stringstream PdbChain::extractSingleResidueFromRecordSection(std::stringstream &pdbFileStream, std::string line)
{
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::stringstream singleResidueSection;
    pdb::ResidueId residueId(line);
    pdb::ResidueId initialResidueId = residueId;
    while(residueId == initialResidueId)
    {
        singleResidueSection << line << std::endl;
        previousLinePosition = pdbFileStream.tellg(); // Save current line position.
        if(!std::getline(pdbFileStream, line))
        {
            break; // // If we hit the end, time to leave.
        }
        residueId = ResidueId(line);
    }
    pdbFileStream.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    gmml::log(__LINE__,__FILE__,gmml::INF, "Single residue section is:\n" + singleResidueSection.str() + "\nEnd of single residue section.");
    return singleResidueSection;
}

void PdbChain::InsertCap(const PdbResidue& refResidue, const std::string& type)
{
    // This approach is bad, should be using templates. When parameter manager is good we can use that to remove the get_carestian stuff
    using GeometryTopology::Coordinate;
    if (type == "NHCH3") // NME
    {
//        int sequenceNumber = refResidue.GetSequenceNumber() + 1; // Single gaps will end up with the same ACE NME resid numbers. Otherwise good.
        const Coordinate* cCoordProtein = refResidue.FindAtom("C")->getCoordinate();
        const Coordinate* caCoordProtein = refResidue.FindAtom("CA")->getCoordinate();
        const Coordinate* oCoordProtein = refResidue.FindAtom("O")->getCoordinate();
        Coordinate nCoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordProtein, caCoordProtein, cCoordProtein, 120.0, 180.0, 1.4);
        Coordinate hCoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordProtein, caCoordProtein, &nCoordNME, 109.0, 180.0, 1.0);
        Coordinate ch3CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, cCoordProtein, &nCoordNME, 125.0, 180.0, 1.48);
        Coordinate hh31CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 180.0, 1.09);
        Coordinate hh32CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 60.0, 1.09);
        Coordinate hh33CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, -60.0, 1.09);
        cds::Residue* newNMEResidue = this->insertNewResidue(std::make_unique<PdbResidue>("NME", &refResidue), refResidue );
        newNMEResidue->addAtom(std::make_unique<PdbAtom>("N", nCoordNME));
        newNMEResidue->addAtom(std::make_unique<PdbAtom>("H", hCoordNME));
        newNMEResidue->addAtom(std::make_unique<PdbAtom>("CH3", ch3CoordNME));
        newNMEResidue->addAtom(std::make_unique<PdbAtom>("HH31", hh31CoordNME));
        newNMEResidue->addAtom(std::make_unique<PdbAtom>("HH32", hh32CoordNME));
        newNMEResidue->addAtom(std::make_unique<PdbAtom>("HH33", hh33CoordNME));
        static_cast<PdbResidue*>(newNMEResidue)->AddTerCard();
    }
    else if (type == "COCH3") // ACE
    {
//        int sequenceNumber = refResidue.GetSequenceNumber() - 1; // Single gaps will end up with the same ACE NME resid numbers. Otherwise good.
        // These are the atoms in residue that I use to build the ACE out from.
        const Coordinate* cCoordProtein = refResidue.FindAtom("C")->getCoordinate();
        const Coordinate* caCoordProtein = refResidue.FindAtom("CA")->getCoordinate();
        const Coordinate* nCoordProtein = refResidue.FindAtom("N")->getCoordinate();
        // This is bad, should use templates loaded from lib/prep file instead.
        Coordinate cCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(cCoordProtein, caCoordProtein, nCoordProtein, 120.0, -130.0, 1.4);
        Coordinate oCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, &cCoordACE, 120.0, 0.0, 1.23);
        Coordinate ch3CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, &cCoordACE, 125.0, 180.0, 1.48);
        Coordinate hh31CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 180.0, 1.09);
        Coordinate hh32CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 60.0, 1.09);
        Coordinate hh33CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, -60.0, 1.09);
        // Ok this next bit is convoluted, but I look up the position of the first atom in the protein residue and insert the new Atom before it, and get passed back the position of the newly created atom, so I can use that when creating the next one and so on.
        // With ACE we want to insert before the residue, so I'm finding the residue before here:
        auto refPosition = this->findPositionOfResidue(&refResidue);
        --refPosition;
        PdbResidue* previousResidue = static_cast<PdbResidue*>((*refPosition).get()); // Its an iterator to a unique ptr, so deref and get the raw. Ugh.
        cds::Residue* newACEResidue = this->insertNewResidue(std::make_unique<PdbResidue>("ACE", previousResidue), *previousResidue);
        newACEResidue->addAtom(std::make_unique<PdbAtom>("C", cCoordACE));
        newACEResidue->addAtom(std::make_unique<PdbAtom>("O", oCoordACE));
        newACEResidue->addAtom(std::make_unique<PdbAtom>("CH3", ch3CoordACE));
        newACEResidue->addAtom(std::make_unique<PdbAtom>("HH31", hh31CoordACE));
        newACEResidue->addAtom(std::make_unique<PdbAtom>("HH32", hh32CoordACE));
        newACEResidue->addAtom(std::make_unique<PdbAtom>("HH33", hh33CoordACE));
        gmml::log(__LINE__, __FILE__, gmml::INF, "Created ACE residue: " + static_cast<PdbResidue*>(newACEResidue)->printId());
    }
}

bool PdbChain::ModifyTerminal(const std::string& type)
{
    if (type == "NH3+") // For now, leaving it to tleap to add the correct H's
    {
        PdbResidue* nTermResidue = this->getNTerminal();
        if (nTermResidue == nullptr) { return false; }
        gmml::log(__LINE__,__FILE__,gmml::INF, "Modifying N Terminal of : " + nTermResidue->printId());
        const PdbAtom* atom = static_cast<const PdbAtom*>(nTermResidue->FindAtom("H"));
        if (atom != nullptr)
        {
            gmml::log(__LINE__,__FILE__,gmml::INF, "Deleting atom with id: " + atom->GetId());
            nTermResidue->deleteAtom(atom);
        }
    }
    else if (type == "CO2-")
    {
        PdbResidue* cTermResidue = this->getCTerminal();
        if (cTermResidue == nullptr) { return false; }
        gmml::log(__LINE__,__FILE__,gmml::INF, "Modifying C Terminal of : " + cTermResidue->printId());
        const cds::Atom* atom = cTermResidue->FindAtom("OXT");
        if (atom == nullptr)
        {
            // I don't like this, but at least it's somewhat contained:
            const cds::Atom* atomCA = cTermResidue->FindAtom("CA");
            const cds::Atom* atomC = cTermResidue->FindAtom("C");
            const cds::Atom* atomO = cTermResidue->FindAtom("O");
            GeometryTopology::Coordinate oxtCoord = GeometryTopology::get_cartesian_point_from_internal_coords(atomCA->getCoordinate(), atomC->getCoordinate(), atomO->getCoordinate(), 120.0, 180.0, 1.25);
            cTermResidue->addAtom(std::make_unique<PdbAtom>("OXT", oxtCoord));
            gmml::log(__LINE__,__FILE__,gmml::INF, "Created new atom named OXT after " + static_cast<const PdbAtom*>(atomO)->GetId());
        }
        else
        {
            gmml::log(__LINE__,__FILE__,gmml::INF, "OXT atom already exists: " + static_cast<const PdbAtom*>(atom)->GetId());
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
        return false;
    }
    return true;
}

// Only makes sense for proteins.
// Assumes vector is populated from N-terminal to C-terminal.
pdb::PdbResidue* PdbChain::getNTerminal()
{
    std::vector<cds::Residue*> proteinResidues = this->getResidues(biology::proteinResidueNames);
    if (proteinResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Looked for terminal residue of chain with protein residues.");
        return nullptr;
    }
    return static_cast<PdbResidue*>(proteinResidues.front());
}

pdb::PdbResidue* PdbChain::getCTerminal()
{
    std::vector<cds::Residue*> proteinResidues = this->getResidues(biology::proteinResidueNames);
    if (proteinResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Looked for terminal residue of chain with protein residues.");
        return nullptr;
    }
    return static_cast<PdbResidue*>(proteinResidues.back());
}


//void PdbChain::addCapsToGaps(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions)
//{
////    // Missing Residues (gaps)
////    gmml::log(__LINE__, __FILE__, gmml::INF, "Gaps");
//// //   std::string previousChainId = "AUniqueInitialString";
////    int previousSequenceNumber = -999999;
//// //   int previousModelNumber = -999999;
////    pdb::PdbResidue* previous = nullptr;
////    for(auto &residue : this->getResidues())
////    {
////        // WE WILL HAVE TO CHECK DISTANCES!!! 1UCY has reverse ordered insertion codes
////        // KABAT can mean skipped numbers that are bonded.
////        if ((previousSequenceNumber != (residue->getNumber() - 1)) && (true) )
////        {
////            gmml::log(__LINE__, __FILE__, gmml::INF, inputOptions.gapNTermination_ + " cap for : " + previous->getId());
////            gmml::log(__LINE__, __FILE__, gmml::INF, inputOptions.gapCTermination_ + " cap for : " + residue->getId());
////            this->InsertCap(*previous, inputOptions.gapCTermination_);
////            this->InsertCap(*residue, inputOptions.gapNTermination_);
////            std::stringstream residueBefore, residueAfter;
////            residueBefore << previous->getNumber() << previous->getInsertionCode();
////            residueAfter << residue->getNumber() << residue->getInsertionCode();
////            ppInfo.missingResidues_.emplace_back(residue->getChainId(), residueBefore.str(), residueAfter.str(), inputOptions.gapCTermination_, inputOptions.gapNTermination_);
////        }
////        previous = residue;
////        previousSequenceNumber = residue->getNumber();
////    }
//}


////////////////////////////////////////////////////////////
////                          MUTATOR                     //
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////////
//void PdbChain::Print(std::ostream &out) const
//{
//    out << "Printing chain " << this->GetChainId() << "\n";
//    for(auto &residue : pdbResidues_)
//    {
//        residue->Print();
//    }
//}
void PdbChain::Write(std::ostream& stream) const
{
    this->WritePdb(stream);
    return;
}
