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
#include "includes/InputSet/PdbFile/pdbModel.hpp"
#include "includes/InputSet/PdbFile/pdbChain.hpp"
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"

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
        std::vector<std::string> coordSectionCards {"MODEL", "ATOM", "ANISOU", "TER", "HETATM", "CONECT"};
        std::vector<std::string> databaseCards {"DBREF", "DBREF1", "DBREF2"};
        if(std::find(coordSectionCards.begin(), coordSectionCards.end(), recordName) != coordSectionCards.end())
        {
            std::stringstream recordSection = this->ExtractHeterogenousRecordSection(pdbFileStream, line, coordSectionCards);
            PdbModel temp = PdbModel(recordSection);
            this->addModel(temp);
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
        else if( (std::find(databaseCards.begin(), databaseCards.end(), recordName) != databaseCards.end()) )
        {
            std::stringstream databaseSection = this->ExtractHeterogenousRecordSection(pdbFileStream, line, databaseCards);
            while(getline(databaseSection, line))
            {
                databaseReferences_.emplace_back(line);
            }
        }
//        else if(recordName == "CONECT")
//        {
//            std::stringstream conectSection = this->ExtractHomogenousRecordSection(pdbFileStream, line, recordName);
//            std::string line;
//            while(getline(conectSection, line))
//            {
//                conectRecords_.emplace_back(line, coordinateSection_); // this is special as it needs the coordinateSection.
//            }
//        }
    }
    return;
}
// Initializers used by constructors
// Should extract all lines that start with the strings in recordNames.
// Returns when it hits a line that does not start with one of those records.
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

//void PdbFile::AddConnection(AtomRecord* atom1, AtomRecord* atom2)
//{
//    conectRecords_.emplace_back(std::vector{atom1, atom2});
//    return;
//}

//void PdbFile::DeleteAtomRecord(AtomRecord* atom)
//{
//    this->GetCoordinateSection().DeleteAtomRecord(atom);
//    return;
//}

//void PdbFile::ModifyNTerminal(const std::string& type, PdbResidue* residue)
//{
//    gmml::log(__LINE__,__FILE__,gmml::INF, "Modifying N Terminal of : " + residue->GetId());
//    if (type == "NH3+")
//    {
//        AtomRecord* atom = residue->FindAtom("H");
//        if (atom != nullptr)
//        {
//            gmml::log(__LINE__,__FILE__,gmml::INF, "Deleting atom with id: " + atom->GetId());
//            residue->DeleteAtomRecord(atom);
//        }
//    }
//    else
//    {
//        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
//    }
//    return;
//}
//
//void PdbFile::ModifyCTerminal(const std::string& type, PdbResidue* residue)
//{
//    gmml::log(__LINE__,__FILE__,gmml::INF, "Modifying C Terminal of : " + residue->GetId());
//    if (type == "CO2-")
//    {
//        AtomRecord* atom = residue->FindAtom("OXT");
//        if (atom == nullptr)
//        {
//            // I don't like this, but at least it's somewhat contained:
//            AtomRecord* atomCA = residue->FindAtom("CA");
//            AtomRecord* atomC = residue->FindAtom("C");
//            AtomRecord* atomO = residue->FindAtom("O");
//            GeometryTopology::Coordinate oxtCoord = GeometryTopology::get_cartesian_point_from_internal_coords(atomCA->GetCoordinate(), atomC->GetCoordinate(), atomO->GetCoordinate(), 120.0, 180.0, 1.25);
//            residue->createAtom("OXT", oxtCoord);
//            gmml::log(__LINE__,__FILE__,gmml::INF, "Created new atom named OXT after " + atomO->GetId());
//        }
//    }
//    else
//    {
//        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
//    }
//    return;
//}

//void PdbFile::SerializeAtomRecordNumbers()
//{
//    for(auto &atomRecord : this->GetCoordinateSection().GetA)
//}

//void PdbFile::InsertCap(const PdbResidue& refResidue, const std::string& type)
//{
//    // This approach is bad, should be using templates. When parameter manager is good we can use that to remove the get_carestian stuff
//    using GeometryTopology::Coordinate;
//    if (type == "NHCH3") // NME
//    {
////        int sequenceNumber = refResidue.GetSequenceNumber() + 1; // Single gaps will end up with the same ACE NME resid numbers. Otherwise good.
//        const Coordinate& cCoordProtein = refResidue.FindAtom("C")->GetCoordinate();
//        const Coordinate& caCoordProtein = refResidue.FindAtom("CA")->GetCoordinate();
//        const Coordinate& oCoordProtein = refResidue.FindAtom("O")->GetCoordinate();
//        Coordinate nCoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordProtein, caCoordProtein, cCoordProtein, 120.0, 180.0, 1.4);
//        Coordinate hCoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordProtein, caCoordProtein, nCoordNME, 109.0, 180.0, 1.0);
//        Coordinate ch3CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, cCoordProtein, nCoordNME, 125.0, 180.0, 1.48);
//        Coordinate hh31CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 180.0, 1.09);
//        Coordinate hh32CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 60.0, 1.09);
//        Coordinate hh33CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, -60.0, 1.09);
//        //AtomRecordIterator atomPosition = this->GetCoordinateSection().FindPositionOfAtom(refResidue.GetLastAtom());
//        PdbResidue *newNMEResidue = this->GetCoordinateSection().CreateNewResidue("NME", "N", nCoordNME, refResidue);
//        newNMEResidue->CreateAtom("H", hCoordNME);
//        newNMEResidue->CreateAtom("CH3", ch3CoordNME);
//        newNMEResidue->CreateAtom("HH31", hh31CoordNME);
//        newNMEResidue->CreateAtom("HH32", hh32CoordNME);
//        newNMEResidue->CreateAtom("HH33", hh33CoordNME);
//        newNMEResidue->AddTerCard();
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("N", "NME", sequenceNumber, nCoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("H", "NME", sequenceNumber, hCoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("CH3", "NME", sequenceNumber, ch3CoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH31", "NME", sequenceNumber, hh31CoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH32", "NME", sequenceNumber, hh32CoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH33", "NME", sequenceNumber, hh33CoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//    }
//    else if (type == "COCH3") // ACE
//    {
////        int sequenceNumber = refResidue.GetSequenceNumber() - 1; // Single gaps will end up with the same ACE NME resid numbers. Otherwise good.
//        // These are the atoms in residue that I use to build the ACE out from.
//        const Coordinate& cCoordProtein = refResidue.FindAtom("C")->GetCoordinate();
//        const Coordinate& caCoordProtein = refResidue.FindAtom("CA")->GetCoordinate();
//        const Coordinate& nCoordProtein = refResidue.FindAtom("N")->GetCoordinate();
//        // This is bad, should use templates loaded from lib/prep file instead.
//        Coordinate cCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(cCoordProtein, caCoordProtein, nCoordProtein, 120.0, -130.0, 1.4);
//        Coordinate oCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, cCoordACE, 120.0, 0.0, 1.23);
//        Coordinate ch3CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, cCoordACE, 125.0, 180.0, 1.48);
//        Coordinate hh31CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 180.0, 1.09);
//        Coordinate hh32CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 60.0, 1.09);
//        Coordinate hh33CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, -60.0, 1.09);
//        // Ok this next bit is convoluted, but I look up the position of the first atom in the protein residue and insert the new Atom before it, and get passed back the position of the newly created atom, so I can use that when creating the next one and so on.
//      //  AtomRecordIterator atomPosition = this->GetCoordinateSection().FindPositionOfAtom(refResidue.GetFirstAtom());
//
//        // With ACE we want to insert before the residue, so I'm finding the residue before here:
//        auto refPosition = this->GetCoordinateSection().FindPositionOfResidue(&refResidue);
//        --refPosition;
//        PdbResidue* previousResidue = (*refPosition).get(); // Its an iterator to a unique ptr, so deref and get the raw. Ugh.
//        PdbResidue *newACEResidue = this->GetCoordinateSection().CreateNewResidue("ACE", "C", cCoordACE, *previousResidue);
//        newACEResidue->CreateAtom("O", oCoordACE);
//        newACEResidue->CreateAtom("CH3", ch3CoordACE);
//        newACEResidue->CreateAtom("HH31", hh31CoordACE);
//        newACEResidue->CreateAtom("HH32", hh32CoordACE);
//        newACEResidue->CreateAtom("HH33", hh33CoordACE);
//
////        --atomPosition; // Want to insert before the first atom, inserting at begin() position is fine.
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("C", "ACE", sequenceNumber, cCoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("O", "ACE", sequenceNumber, oCoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("CH3", "ACE", sequenceNumber, ch3CoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH31", "ACE", sequenceNumber, hh31CoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH32", "ACE", sequenceNumber, hh32CoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
////        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH33", "ACE", sequenceNumber, hh33CoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        gmml::log(__LINE__, __FILE__, gmml::INF, "Created ACE residue: " + newACEResidue->GetId());
//    }
//}

pdb::PreprocessorInformation PdbFile::PreProcess(PreprocessorOptions inputOptions)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Preprocesssing has begun");
    PreprocessorInformation ppInfo;
    // CYS Disulfide bonds
    gmml::log(__LINE__, __FILE__, gmml::INF, "Cys disulphide bonds");
    for (auto &models: this->GetCoordinateSection().GetModels())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "I have the models");
        std::vector<pdb::PdbResidue*> cysResidues;
        for(auto &residue : models)
        {
            if (residue->GetName() == "CYS" || residue->GetName() == "CYX")
            {
                cysResidues.push_back(residue);
            }
        }
        if (cysResidues.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "No CYS or CYX residues detected in this structure\n");
        }
        for (std::vector<pdb::PdbResidue*>::iterator it1 = cysResidues.begin(); it1 != cysResidues.end(); ++it1)
        {
            PdbResidue* res1 = *it1;
            AtomRecord* sg1 = res1->FindAtom("SG");
            for (std::vector<pdb::PdbResidue*>::iterator it2 = std::next(it1, 1); it2 != cysResidues.end(); ++it2)
            {
                PdbResidue* res2 = *it2;
                AtomRecord* sg2 = res2->FindAtom("SG");
                if ( (sg1 != nullptr) && (sg2 != nullptr) )
                {
                    double distance = sg1->CalculateDistance(sg2);
                    if (distance < gmml::dSulfurCutoff && distance > 0.001)
                    {
                        res1->SetName("CYX");
                        res2->SetName("CYX");
                        this->AddConnection(res1->FindAtom("SG"), res2->FindAtom("SG"));
                        ppInfo.cysBondResidues_.emplace_back(res1->GetId(), res2->GetId(), distance);
                        std::stringstream message;
                        message << "Bonding " << res1->GetId() << " and " << res2->GetId() << " with distance " << distance;
                        gmml::log(__LINE__, __FILE__, gmml::INF, message.str());
                    }
                }
            }
        }
    }
    // HIS protonation, user specified:
    gmml::log(__LINE__, __FILE__, gmml::INF, "User His protonation");
    for(auto &userSelectionPair : inputOptions.hisSelections_)
    {
        this->GetCoordinateSection().ChangeResidueName(userSelectionPair.first, userSelectionPair.second);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Auto His protonation");
    // HIS protonation, automatic handling.
    for(auto &residue : this->GetCoordinateSection().GetResidues())
    {
        if (residue->GetName() == "HIE" || residue->GetName() == "HID" || residue->GetName() == "HIP")
        {
            ppInfo.hisResidues_.emplace_back(residue->GetId());
        }
        else if (residue->GetName() == "HIS")
        {
            if ( (residue->FindAtom("HE2") == nullptr) && (residue->FindAtom("HD1") != nullptr) )
            {
                residue->SetName("HID");
            }
            else if ( (residue->FindAtom("HE2") != nullptr) && (residue->FindAtom("HD1") != nullptr) )
            {
                residue->SetName("HIP");
            }
            else // HIE is default
            {
                residue->SetName("HIE");
            }
            gmml::log(__LINE__, __FILE__, gmml::INF, "About to emplaceBack Id");
            gmml::log(__LINE__, __FILE__, gmml::INF, residue->GetId());
            ppInfo.hisResidues_.emplace_back(residue->GetId());
        }
    }
    //Chain terminations
    gmml::log(__LINE__, __FILE__, gmml::INF, "Chain terminations");
    for (std::vector<pdb::PdbResidue*> chainOfResidues : this->GetCoordinateSection().GetProteinChains())
    {
       //Do the thing
       this->ModifyNTerminal(inputOptions.chainNTermination_, chainOfResidues.front());
       this->ModifyCTerminal(inputOptions.chainCTermination_, chainOfResidues.back());
       //Log the thing
       gmml::log(__LINE__, __FILE__, gmml::INF, "N term : " + chainOfResidues.front()->GetId());
       gmml::log(__LINE__, __FILE__, gmml::INF, "C term : " + chainOfResidues.back()->GetId());
       //Report the thing
       std::stringstream startIndex, endIndex;
       startIndex << chainOfResidues.front()->GetSequenceNumber() << chainOfResidues.front()->GetInsertionCode();
       endIndex << chainOfResidues.back()->GetSequenceNumber() << chainOfResidues.back()->GetInsertionCode();
       ppInfo.chainTerminals_.emplace_back(chainOfResidues.front()->GetChainId(), startIndex.str(), endIndex.str(), inputOptions.chainNTermination_, inputOptions.chainCTermination_);
    }
    // Missing Residues (gaps)
    gmml::log(__LINE__, __FILE__, gmml::INF, "Gaps");
    std::string previousChainId = "AUniqueInitialString";
    int previousSequenceNumber = -999999;
    int previousModelNumber = -999999;
    pdb::PdbResidue* previous = nullptr;
    for(auto &residue : this->GetCoordinateSection().GetResidues())
    {
        if ((previousSequenceNumber != (residue->GetSequenceNumber() - 1)) && (previousChainId == residue->GetChainId() ) && (previousModelNumber == residue->GetModelNumber()) )
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, inputOptions.gapNTermination_ + " cap for : " + previous->GetId());
            gmml::log(__LINE__, __FILE__, gmml::INF, inputOptions.gapCTermination_ + " cap for : " + residue->GetId());
            this->InsertCap(*previous, inputOptions.gapCTermination_);
            this->InsertCap(*residue, inputOptions.gapNTermination_);
            std::stringstream residueBefore, residueAfter;
            residueBefore << previous->GetSequenceNumber() << previous->GetInsertionCode();
            residueAfter << residue->GetSequenceNumber() << residue->GetInsertionCode();
            ppInfo.missingResidues_.emplace_back(residue->GetChainId(), residueBefore.str(), residueAfter.str(), inputOptions.gapCTermination_, inputOptions.gapNTermination_);
        }
        previous = residue;
        previousChainId = residue->GetChainId();
        previousSequenceNumber = residue->GetSequenceNumber();
        previousModelNumber = residue->GetModelNumber();
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Unrecognized or missing residues and atoms");
    parameters::Manager parmManager;
    for(auto &residue : this->GetCoordinateSection().GetResidues())
    {
        std::vector<std::string> parmAtomNames = parmManager.GetAtomNamesForResidue(residue->GetParmName());
        std::vector<std::string> parmHeavyAtomNames = parmManager.GetHeavyAtomNamesForResidue(residue->GetParmName());
        // Unrecognized residue->
        if (parmAtomNames.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "ParmManager did not recognize residue: " + residue->GetParmName());
            ppInfo.unrecognizedResidues_.emplace_back(residue->GetId());
        }
        else // Recognized residue->
        {
            std::vector<std::string> pdbAtomNames = residue->GetAtomNames();
            for (auto &parmHeavyAtomName : parmHeavyAtomNames) // What heavy atoms are missing from the pdb residue?
            {
                if ( std::find(pdbAtomNames.begin(), pdbAtomNames.end(), parmHeavyAtomName) == pdbAtomNames.end() )
                { // Residue missing a heavy atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Atom named " + parmHeavyAtomName + " missing from " + residue->GetId());
                    ppInfo.missingHeavyAtoms_.emplace_back(parmHeavyAtomName, residue->GetId());
                }
            }
            for (auto &pdbAtomName : pdbAtomNames) // What atoms in the pdb residue are unrecognized?
            {
                if ( std::find(parmAtomNames.begin(), parmAtomNames.end(), pdbAtomName) == parmAtomNames.end() )
                {
                    // Residue contains unrecognized atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Unrecognized atom named " + pdbAtomName + " in " + residue->GetId());
                    ppInfo.unrecognizedAtoms_.emplace_back(pdbAtomName, residue->GetId());
                }
            }
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Preprocessing completed");
    return ppInfo;
}

//void PdbFile::Print(std::ostream& out) const
//{
//    coordinateSection_.Print(out);
//}
void PdbFile::Write(const std::string outName) const
{
    std::ofstream outFileStream;
    try
    {
        outFileStream.open(outName.c_str());
    }
    catch(...)
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Output file could not be created:\n" + outName);
        throw std::runtime_error("Output file could not be created:\n" + outName);
    }
    try
    {
        this->Write(outFileStream);
    }
    catch(...)
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Error when writing pdbFile class to file:\n" + outName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + outName);
    }
}
void PdbFile::Write(std::ostream& out) const
{
    this->GetHeaderRecord().Write(out);
    this->GetTitleRecord().Write(out);
    this->GetAuthorRecord().Write(out);
    this->GetJournalRecord().Write(out);
    this->GetRemarkRecord().Write(out);
    for (auto &dbref : this->GetDatabaseReferences())
    {
        dbref.Write(out);
    }
    for (auto &model : this->getModels())
    {
        model->Write(out);
    }
    for (auto &conect : this->GetConectRecords())
    {
        conect.Write(out);
    }
}
