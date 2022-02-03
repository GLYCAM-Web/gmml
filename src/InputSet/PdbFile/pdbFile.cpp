#include <fstream>      // std::ifstream
#include <algorithm>    // std::find
#include "includes/InputSet/PdbFile/pdbFile.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/InputSet/PdbFile/databaseReferenceRecord.hpp"
#include "includes/common.hpp" // gmml::dSulfurCutoff
#include "includes/ParameterSet/parameterManager.hpp" // for preprocssing
#include "includes/GeometryTopology/geometrytopology.hpp"

using pdb::PdbFile;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFile::PdbFile()
{
    //this->Initialize();
    inFilePath_ = "GMML-Generated";
}
PdbFile::PdbFile(const std::string &pdbFilePath) : inFilePath_(pdbFilePath)
{
    std::ifstream pdbFileStream(pdbFilePath);
    if(pdbFileStream.fail())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Could not open this file: " + pdbFilePath);
        throw std::runtime_error("PdbFile constructor could not open this file: " + pdbFilePath);
    }
    gmml::log(__LINE__, __FILE__,  gmml::INF, "File opened: " + pdbFilePath + ". Ready to parse!");
    this->ParseInFileStream(pdbFileStream);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Finished parsing " + pdbFilePath);
}
void PdbFile::ParseInFileStream(std::ifstream& pdbFileStream)
{
    std::cout << "Parsing inputfile\n";
    for (std::string line; std::getline(pdbFileStream, line); )
    {
        //std::cout << "Parsing the line: " << line << "\n";
        codeUtils::ExpandLine(line, pdb::iPdbLineLength);
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        std::vector<std::string> cards {"MODEL", "ATOM", "ANISOU", "TER", "HETATM", "ENDMDL"};
        if(std::find(cards.begin(), cards.end(), recordName) != cards.end())
        {
            std::stringstream recordSection = this->ExtractHeterogenousRecordSection(pdbFileStream, line, cards);
            coordinateSection_ = CoordinateSection(recordSection);
        }
        else if(recordName == "HEADER")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            headerRecord_ = HeaderRecord(recordSection);
        }
        else if(recordName == "TITLE")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            titleRecord_ = TitleRecord(recordSection);
        }
        else if(recordName == "AUTHOR")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            authorRecord_ = AuthorRecord(recordSection);
        }
        else if(recordName == "JRNL")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            journalRecord_ = JournalRecord(recordSection);
        }
        else if(recordName == "REMARK")
        {
            std::stringstream recordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            remarkRecord_ = RemarkRecord(recordSection);
        }
        else if(recordName == "CONECT")
        {
            std::stringstream conectSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
            std::string line;
            while(getline(conectSection, line))
            {
                conectRecords_.emplace_back(line, coordinateSection_); // this is special as it needs the coordinateSection.
            }
        }
    }
    return;
}
// Initializers used by constructors
std::stringstream PdbFile::ExtractHeterogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, const std::vector<std::string> recordNames)
{
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::stringstream recordSection;
    std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
    while(std::find(recordNames.begin(), recordNames.end(), recordName) != recordNames.end())
    {
        std::stringstream partialRecordSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
        recordSection << partialRecordSection.str();
        previousLinePosition = pdbFileStream.tellg(); // Save current line position.
        if(!std::getline(pdbFileStream, line)) // If we hit the end
        {
            break; // Time to leave.
        }
        recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
    }
    pdbFileStream.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    //std::cout << "Returning this hetero section:\n" << heteroRecordSection.str() << "\nThe End of Hetero Section.\n";
    return recordSection;
}

// Goes through a section of the PDB file that contains the same header section. E.g. HEADER.
// If the header changes, it goes back to the previous line. I wanted the out while loop to trigger the new line. This means I don't have to check recordName between if statement and can have else if.
std::stringstream PdbFile::ExtractHomogenousRecordSection(std::ifstream &pdbFileStream, std::string &line, std::string recordName)
{
    std::stringstream recordSection;
    codeUtils::ExpandLine(line, pdb::iPdbLineLength);
    recordSection << line << std::endl;
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::string previousName = recordName;
    while ( (std::getline(pdbFileStream, line)) )
    {
        codeUtils::ExpandLine(line, pdb::iPdbLineLength);
        recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        if (recordName == previousName)
        {
            recordSection << line << std::endl;
            previousName = recordName;
            previousLinePosition = pdbFileStream.tellg(); // Save current line position.
        }
        else
        {
            break;
        }
    }
    pdbFileStream.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    //std::cout << "At end. Returning this record section:\n" << recordSection.str() << "\nEND RECORD SECTION\n";
    return recordSection;
}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
std::string PdbFile::GetUniprotIDs() const
{
    std::string UniprotIDs = "";
    for (auto &databaseReference : this->GetDatabaseReferences())
    {
        UniprotIDs += databaseReference.GetUniprotID();
    }
    return UniprotIDs;
}

const float& PdbFile::GetResolution() const
{
    return this->GetRemarkRecord().GetResolution();
}

const float& PdbFile::GetBFactor() const
{
    return this->GetRemarkRecord().GetBFactor();
}

void PdbFile::AddConnection(AtomRecord* atom1, AtomRecord* atom2)
{
    conectRecords_.emplace_back(std::vector{atom1, atom2});
    return;
}

void PdbFile::DeleteAtomRecord(AtomRecord* atom)
{
    this->GetCoordinateSection().DeleteAtomRecord(atom);
    return;
}

void PdbFile::ModifyNTerminal(const std::string& type, PdbResidue* residue)
{
    gmml::log(__LINE__,__FILE__,gmml::INF, "Modifying N Terminal of : " + residue->GetId());
    if (type == "Zwitterionic")
    {
        AtomRecord* atom = residue->FindAtom("H");
        if (atom != nullptr)
        {
            gmml::log(__LINE__,__FILE__,gmml::INF, "Deleting atom with id: " + atom->GetId());
            this->DeleteAtomRecord(atom);
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}

void PdbFile::ModifyCTerminal(const std::string& type, PdbResidue* residue)
{
    gmml::log(__LINE__,__FILE__,gmml::INF, "Modifying C Terminal of : " + residue->GetId());
    if (type == "Zwitterionic")
    {
        AtomRecord* atom = residue->FindAtom("OXT");
        if (atom == nullptr)
        {
            // I don't like this, but at least it's somewhat contained:
            AtomRecord* atomCA = residue->FindAtom("CA");
            AtomRecord* atomC = residue->FindAtom("C");
            AtomRecord* atomO = residue->FindAtom("O");
            GeometryTopology::Coordinate oxtCoord = GeometryTopology::get_cartesian_point_from_internal_coords(atomCA->GetCoordinate(), atomC->GetCoordinate(), atomO->GetCoordinate(), 120.0, 180.0, 1.25);
            this->GetCoordinateSection().CreateNewAtomRecord("OXT", oxtCoord, atomO);
            gmml::log(__LINE__,__FILE__,gmml::INF, "Created new atom named OXT after " + atomO->GetId());
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}

void PdbFile::InsertCap(const PdbResidue& residue, const std::string& type)
{
    // This approach is bad, should be using templates. When parameter manager is good we can use that to remove the get_carestian stuff
    using GeometryTopology::Coordinate;
    if (type == "ACE")
    {
        int sequenceNumber = residue.GetSequenceNumber() - 1; // Single gaps will end up with the same ACE NME resid numbers. Otherwise good.

        const Coordinate& cCoordProtein = residue.FindAtom("C")->GetCoordinate();
        const Coordinate& caCoordProtein = residue.FindAtom("CA")->GetCoordinate();
        const Coordinate& nCoordProtein = residue.FindAtom("N")->GetCoordinate();

        Coordinate cCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(cCoordProtein, caCoordProtein, nCoordProtein, 120.0, -130.0, 1.4);
        Coordinate oCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, cCoordACE, 120.0, 0.0, 1.23);
        Coordinate ch3CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, cCoordACE, 125.0, 180.0, 1.48);
        Coordinate hh31CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 180.0, 1.09);
        Coordinate hh32CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 60.0, 1.09);
        Coordinate hh33CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, -60.0, 1.09);
        AtomRecord* cAtomACE = this->GetCoordinateSection().CreateNewAtomRecord("C", "ACE", sequenceNumber, cCoordACE, residue.GetChainId(), residue.GetModelNumber(), residue.GetLastAtom());
        AtomRecord* oAtomACE = this->GetCoordinateSection().CreateNewAtomRecord("O", "ACE", sequenceNumber, oCoordACE, residue.GetChainId(), residue.GetModelNumber(), cAtomACE);
        AtomRecord* ch3AtomACE = this->GetCoordinateSection().CreateNewAtomRecord("CH3", "ACE", sequenceNumber, ch3CoordACE, residue.GetChainId(), residue.GetModelNumber(), oAtomACE);
//        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
//        {
//            PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
//            int index = distance(ordered_atoms_of_residue.begin(), it3);
//            PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number,
//                    atom_of_residue->GetAtomName(),
//                    atom_of_residue->GetAtomAlternateLocation(),
//                    atom_of_residue->GetAtomResidueName(),
//                    atom_of_residue->GetAtomChainId(),
//                    sequence_number,
//                    atom_of_residue->GetAtomInsertionCode(),
//                    coordinate_set.at(index),
//                    atom_of_residue->GetAtomOccupancy(),
//                    atom_of_residue->GetAtomTempretureFactor(),
//                    atom_of_residue->GetAtomElementSymbol(),
//                    atom_of_residue->GetAtomCharge(),
//                    atom->GetAlternateAtomCards());
//            updated_atoms[serial_number] = new_atom;
//            updated_atoms_vector.push_back(new_atom);
//            serial_number++;
//        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Created ACE residue numbered: " + std::to_string(sequenceNumber));
    }
}

pdb::PreprocessorInformation PdbFile::PreProcess(PreprocessorOptions options)
{
    std::cout << "Preprocesssing has begun\n";
    PreprocessorInformation ppInfo;
    // CYS Disulfide bonds
    std::vector<pdb::PdbResidue> cysResidues = this->GetCoordinateSection().FindResidues("CYS");
    std::vector<pdb::PdbResidue> cyxResidues = this->GetCoordinateSection().FindResidues("CYX");
    cysResidues.insert( cysResidues.end(), cyxResidues.begin(), cyxResidues.end()); // concatenates the vectors.
    if (cysResidues.empty())
    {
        std::cout << "No CYS or CYX residue detected in this structure\n";
    }
    for (auto &cysRes1 : cysResidues)
    {
        AtomRecord* sg1 = cysRes1.FindAtom("SG");
        for (auto &cysRes2 : cysResidues)
        {
            AtomRecord* sg2 = cysRes2.FindAtom("SG");
            if ( (sg1 != nullptr) && (sg2 != nullptr) )
            {
                double distance = sg1->CalculateDistance(sg2);
                if (distance < gmml::dSulfurCutoff && distance > 0.001)
                {
                    cysRes1.SetName("CYX");
                    cysRes2.SetName("CYX");
                    this->AddConnection(cysRes1.FindAtom("SG"), cysRes2.FindAtom("SG"));
                    std::cout << "Bonding " << cysRes1.GetId() << " and " << cysRes2.GetId() << " with distance " << distance << "\n";
                }
            }
        }
    }
    // HIS Protonation
    for (auto &hisRes : this->GetCoordinateSection().FindResidues("HIS"))
    {
        // HID residue
        if ( (hisRes.FindAtom("HE2") == nullptr) && (hisRes.FindAtom("HD1") != nullptr) )
        {
            hisRes.SetName("HID");
        }
        // HIP residue
        else if ( (hisRes.FindAtom("HE2") != nullptr) && (hisRes.FindAtom("HD1") != nullptr) )
        {
            hisRes.SetName("HIP");
        }
        // HIE is default
        else
        {
            hisRes.SetName("HIE");
        }
    }
    //Chain terminations
    for (std::vector<pdb::PdbResidue> chainOfResidues : this->GetCoordinateSection().GetProteinChains())
    {
       gmml::log(__LINE__, __FILE__, gmml::INF, "N term : " + chainOfResidues.front().GetId());
       gmml::log(__LINE__, __FILE__, gmml::INF, "C term : " + chainOfResidues.back().GetId());
       this->ModifyNTerminal("Zwitterionic", &chainOfResidues.front()); // Check input for Zwitterionic or something else
       this->ModifyCTerminal("Zwitterionic", &chainOfResidues.back());  // Check input for Zwitterionic or something else
    }
    // Missing Residues (gaps)
    std::string previousChainId = "AUniqueInitialString";
    int previousSequenceNumber = 999999;
    pdb::PdbResidue* previous = nullptr;
    for(auto &residue : this->GetCoordinateSection().GetResidues())
    {
        if ( (previousChainId == residue.GetChainId() ) && (previousSequenceNumber != (residue.GetSequenceNumber() - 1)) )
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "NME cap for : " + previous->GetId());
            gmml::log(__LINE__, __FILE__, gmml::INF, "ACE cap for : " + residue.GetId());
            this->InsertCap(residue, "NME");
            this->InsertCap(residue, "ACE");
        }
        previous = &residue;
        previousChainId = residue.GetChainId();
        previousSequenceNumber = residue.GetSequenceNumber();
    }
    parameters::Manager parmManager;
    for(auto &residue : this->GetCoordinateSection().GetResidues())
    {
        std::vector<std::string> parmAtomNames = parmManager.GetAtomNamesForResidue(residue.GetParmName());
        std::vector<std::string> parmHeavyAtomNames = parmManager.GetHeavyAtomNamesForResidue(residue.GetParmName());
        // Unrecognized residue.
        if (parmAtomNames.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "ParmManager did not recognize residue: " + residue.GetParmName());
        }
        else // Recognized residue.
        {
            std::vector<std::string> pdbAtomNames = residue.GetAtomNames();
            // What heavy atoms are missing from the pdb residue?
            for (auto &parmHeavyAtomName : parmHeavyAtomNames)
            {
                if ( std::find(pdbAtomNames.begin(), pdbAtomNames.end(), parmHeavyAtomName) == pdbAtomNames.end() )
                {
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Atom named " + parmHeavyAtomName + " missing from " + residue.GetId());
                    // Residue missing a heavy atom.
                }
            }
            // What atoms in the pdb residue are unrecognized?
            for (auto &pdbAtomName : pdbAtomNames)
            {
                if ( std::find(parmAtomNames.begin(), parmAtomNames.end(), pdbAtomName) == parmAtomNames.end() )
                {
                    // Residue contains unrecognized atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Unrecognized atom named " + pdbAtomName + " in " + residue.GetId());
                }
            }
        }
    }
    return ppInfo;
}

void PdbFile::Print(std::ostream& out) const
{
    coordinateSection_.Print(out);
}
