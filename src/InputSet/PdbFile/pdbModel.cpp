#include <algorithm> // std::find

#include "includes/InputSet/PdbFile/pdbModel.hpp"
#include "includes/InputSet/PdbFile/pdbChain.hpp"
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/common.hpp" // gmml::PROTEINS
#include "includes/ParameterSet/parameterManager.hpp" // for preprocssing


using pdb::PdbModel;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModel::PdbModel() {}

PdbModel::PdbModel(std::stringstream &stream_block)
{
    int currentModelNumber = 1;
    std::string previousResidueId = "InitialValue";
    std::string line;
    while(getline(stream_block, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        // MODEL cards
        if(recordName == "MODEL")
        {
            try
            {
                currentModelNumber = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(10,4)));
            }
            catch (...)
            {
                gmml::log(__LINE__, __FILE__, gmml::WAR, "Model issue: this ain't an int: " + codeUtils::RemoveWhiteSpace(line.substr(10,4)));
                currentModelNumber = 1; // Seems like a reasonable default.
            }
            this->setNumber(currentModelNumber);
        }
        // ATOM
        else if ( (recordName == "ATOM") || (recordName == "HETATM") )
        {
            // Gimme everything with the same chain, can be everything with no chain.
            // Function that will read from stringstream until chain ID changes or TER or just not ATOM/HETATM
            std::stringstream singleChainSection = this->extractSingleChainFromRecordSection(stream_block, line, this->extractChainId(line));
            //this->addMolecule(PdbChain(singleChainSection, this->extractChainId(line)));
//            gmml::log(__LINE__,__FILE__,gmml::INF, "I'm about to die now?");
            this->addMolecule(std::make_unique<PdbChain>(singleChainSection, this->extractChainId(line)));
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::INF, "PdbModel Constructor Complete Captain");
    return;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

// This is bad as it repeats how to read a line and how to create a residue ID, but I need to know which residue to put the atom into before I construct the atom.
//std::string PdbModel::PeekAtResidueId(const std::string &line)
//{
//    // Dealing with number overruns for serialNumber and residueNumber
//    int shift = codeUtils::GetSizeOfIntInString(line.substr(12));
//    std::string residueName = codeUtils::RemoveWhiteSpace(line.substr(17 + shift, 3));
//    std::string chainId = codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
//    int secondShift = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
//    std::string residueNumber = codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift));
//    // Insertion code gets shifted right by every overrun in residue number.
//    std::string insertionCode = codeUtils::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
//    return residueName + "_" + residueNumber + "_" + insertionCode + "_" + chainId;
//}

std::string PdbModel::extractChainId(const std::string &line)
{   // serialNumber can overrun into position 12 in input.
    int shift = codeUtils::GetSizeOfIntInString(line.substr(12));
    return codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
}

std::stringstream PdbModel::extractSingleChainFromRecordSection(std::stringstream &stream_block, std::string line, const std::string& initialChainID)
{
    std::streampos previousLinePosition = stream_block.tellg(); // Save current line position
    std::stringstream singleChainSection;
    std::string chainID = initialChainID;
    std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
    while((chainID == initialChainID) && (recordName != "TER") )
    {
        singleChainSection << line << std::endl;
        previousLinePosition = stream_block.tellg(); // Save current line position.
        if(!std::getline(stream_block, line))
        {
            break; // // If we hit the end, time to leave.
        }
        chainID = this->extractChainId(line);
        recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
    }
    stream_block.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    gmml::log(__LINE__,__FILE__,gmml::INF, "Single chain section is:\n" + singleChainSection.str() + "\nEnd of single chain section.");
    return singleChainSection;
}

void PdbModel::addConectRecord(const pdbAtom* atom1, const pdbAtom* atom2)
{
    conectRecords_.emplace_back(std::vector{atom1, atom2});
    return;
}

//std::vector<std::vector<pdb::PdbResidue*>> PdbAssembly::GetModels() const
//{
//    std::vector<std::vector<pdb::PdbResidue*>> models;
//    int previousModel = -123456789; // It would need to be very wrong to be this value.
//    for(auto &residue : residues_)
//    {
//        if (residue->GetModelNumber() != previousModel )
//        {
//            models.push_back(std::vector<PdbResidue*>{residue.get()});
//            previousModel = residue->GetModelNumber();
//        }
//        else
//        {
//            models.back().push_back(residue.get());
//        }
//    }
//    return models;
//}

//std::vector<pdb::PdbChain*> PdbModel::GetProteinChains()
//{
//    std::vector<pdb::PdbChain*> proteinChains;
//    std::string previousChain = "InitialValue";
//    for(auto &chain : this->getMolecules())
//    {
//        proteinChains.emplace_back(chain->getResidues(gmml::proteinResidueNames));
//        // This gmml::PROTEINS seems weird, but whatever, it works.
//        if( std::find( gmml::PROTEINS, ( gmml::PROTEINS + gmml::PROTEINSSIZE ), residue->getName() ) != ( gmml::PROTEINS + gmml::PROTEINSSIZE ) )
//        {
//            if ( residue->getChainId() != previousChain || residue->getModelNumber() != previousModelNumber )
//            {
//                proteinChains.emplace_back(std::vector<PdbResidue*>{residue.get()});
//                previousChain = residue->GetChainId();
//                previousModelNumber = residue->GetModelNumber();
//            }
//            else
//            {
//                proteinChains.back().push_back(residue.get());
//            }
//        }
//    }
//    // Labels
//    for(auto &chain : proteinChains)
//    {
//        chain.front()->addLabel("NTerminal");
//        chain.back()->addLabel("CTerminal");
//    }
//    return proteinChains;
//}

//std::vector<pdb::PdbResidue*> PdbAssembly::GetResidues() const
//{
//    std::vector<pdb::PdbResidue*> residues; // reserve the right amount of memory.
//    for(auto &resUniquePtr : residues_)
//    {
//        residues.push_back(resUniquePtr.get());
//    }
//    return residues;
//}

//std::vector<pdb::PdbResidue*> PdbModel::FindResidues(const std::string selector)
//{
//    std::vector<pdb::PdbResidue*> matchingResidues;
//    for(auto &residue : residues_)
//    {
//        std::size_t found = residue->GetId().find(selector);
//        if(found != std::string::npos)
//        {
//            matchingResidues.push_back(residue.get());
//        }
//    }
//    return matchingResidues;
//}
//
void PdbModel::ChangeResidueName(const std::string& selector, const std::string& newName)
{
    for(auto &residue : this->getResidues())
    {
        std::size_t found = residue->printId().find(selector);
        if(found != std::string::npos)
        {
            residue->setName(newName);
            return;
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR, "Could not find residue to rename with this selector " + selector);
    return;
}

const pdb::pdbAtom* PdbModel::FindAtom(const int& serialNumber) const
{
    for(auto &atom : this->getAtoms())
    {
        if (atom->GetSerialNumber() == serialNumber)
        {
            return atom;
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR, "Could not find atom with this serialNumber " + serialNumber);
    return nullptr;
}

//{
//    for(auto &atom : this->getAtoms())
//    {
//        if (atom->GetSerialNumber() == serialNumber)
//        {
//            return atom;
//        }
//    }
//    gmml::log(__LINE__, __FILE__, gmml::WAR, "Could not find atom with this serialNumber " + serialNumber);
//    return nullptr;
//}

//////////////////////////////////////////////////////////
//                      FUNCTIONS                       //
//////////////////////////////////////////////////////////
//void CoordinateSection::SerializeAtomRecordSerialNumbers()
//{
//
//}
void PdbModel::preProcessCysResidues(pdb::PreprocessorInformation &ppInfo)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Start CYS preprocessing for this Model\n");
    std::vector<pdb::PdbResidue*> cysResidues = this->getResiduesWithName(std::vector<std::string> {"CYS", "CYX"});
    if (cysResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "No CYS or CYX residues detected in this structure\n");
    }
    for (std::vector<pdb::PdbResidue*>::iterator it1 = cysResidues.begin(); it1 != cysResidues.end(); ++it1)
    { // I want to go through the list and compare from current item to end. Thus it2 = std::next it1

        PdbResidue* cysRes1 = *it1;
        pdbAtom* sgAtom1 = cysRes1->FindAtom("SG");
        for (std::vector<pdb::PdbResidue*>::iterator it2 = std::next(it1, 1); it2 != cysResidues.end(); ++it2)
        {
            PdbResidue* cysRes2 = *it2;
            pdbAtom* sgAtom2 = cysRes2->FindAtom("SG");
            if ( (sgAtom1 != nullptr) && (sgAtom2 != nullptr) )
            {
                //gmml::log(__LINE__, __FILE__, gmml::INF, "Found SG ATOMS");
                double distance = sgAtom1->calculateDistance(sgAtom2);
                if (distance < gmml::dSulfurCutoff && distance > 0.001)
                {
                    //gmml::log(__LINE__, __FILE__, gmml::INF, "Distance less than cutoff");
                    cysRes1->setName("CYX");
                    cysRes2->setName("CYX");
                    //gmml::log(__LINE__, __FILE__, gmml::INF, "Names set");
                    this->addConectRecord(sgAtom1, sgAtom2);
                    ppInfo.cysBondResidues_.emplace_back(cysRes1->getId(), cysRes2->getId(), distance);
                    //gmml::log(__LINE__, __FILE__, gmml::INF, "ThisNoHappen?");
                    std::stringstream message;
                    message << "Bonding " << cysRes1->printId() << " and " << cysRes2->printId() << " with distance " << distance;
                    gmml::log(__LINE__, __FILE__, gmml::INF, message.str());
                }
            }
        }
    }
    return;
}

void PdbModel::preProcessHisResidues(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions)
{
    // HIS protonation, user specified:
    gmml::log(__LINE__, __FILE__, gmml::INF, "User His protonation");
    for(auto &userSelectionPair : inputOptions.hisSelections_)
    {
        this->ChangeResidueName(userSelectionPair.first, userSelectionPair.second);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Auto His protonation");
    // HIS protonation, automatic handling.
    for(auto &residue : this->getResidues())
    {
        if (residue->getName() == "HIE" || residue->getName() == "HID" || residue->getName() == "HIP")
        {
            ppInfo.hisResidues_.emplace_back(residue->getId());
        }
        else if (residue->getName() == "HIS")
        {
            if ( (residue->FindAtom("HE2") == nullptr) && (residue->FindAtom("HD1") != nullptr) )
            {
                residue->setName("HID");
            }
            else if ( (residue->FindAtom("HE2") != nullptr) && (residue->FindAtom("HD1") != nullptr) )
            {
                residue->setName("HIP");
            }
            else // HIE is default
            {
                residue->setName("HIE");
            }
            gmml::log(__LINE__, __FILE__, gmml::INF, "About to emplaceBack Id");
            gmml::log(__LINE__, __FILE__, gmml::INF, residue->printId());
            ppInfo.hisResidues_.emplace_back(residue->getId());
        }
    }
    return;
}

void PdbModel::preProcessChainTerminals(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Chain terminations");
    for (auto &chain : this->getMolecules())
    {
        gmml::log(__LINE__,__FILE__,gmml::INF, "Chain termination processing started for this chain");
        //Do the thing
        if (chain->ModifyTerminal(inputOptions.chainNTermination_) && chain->ModifyTerminal(inputOptions.chainCTermination_) )
        {
            //Log the thing
            PdbResidue* nTer = chain->getNTerminal();
            PdbResidue* cTer = chain->getCTerminal();
            gmml::log(__LINE__, __FILE__, gmml::INF, "N term : " + nTer->printId());
            gmml::log(__LINE__, __FILE__, gmml::INF, "C term : " + cTer->printId());
            //Report the thing
            ppInfo.chainTerminals_.emplace_back(nTer->getChainId(), nTer->getNumberAndInsertionCode(), cTer->getNumberAndInsertionCode(), inputOptions.chainNTermination_, inputOptions.chainCTermination_);
        }
        else
        {
            gmml::log(__LINE__,__FILE__,gmml::INF, "Could not modify terminals of this chain.");
        }
        gmml::log(__LINE__,__FILE__,gmml::INF, "Preprocessing complete for this chain");
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Chain termination processing complete");
    return;
}

void PdbModel::preProcessGaps(pdb::PreprocessorInformation &ppInfo, const pdb::PreprocessorOptions& inputOptions)
{
    // Missing Residues (gaps)
    gmml::log(__LINE__, __FILE__, gmml::INF, "Gaps");
    std::string previousChainId = "AUniqueInitialString";
    int previousSequenceNumber = -999999;
    //   int previousModelNumber = -999999;
    pdb::PdbResidue* previous = nullptr;
    for(auto &chain : this->getMolecules())
    {
        for(auto &residue : chain->getResidues())
        {   // ToDo WE WILL HAVE TO CHECK DISTANCES!!! 1UCY has reverse ordered insertion codes
            // KABAT can mean skipped numbers that are bonded.
            if ((previousSequenceNumber != (residue->getNumber() - 1)) && previousChainId == residue->getChainId())
            {
                //Log it
                gmml::log(__LINE__, __FILE__, gmml::INF, inputOptions.gapNTermination_ + " cap for : " + previous->printId());
                gmml::log(__LINE__, __FILE__, gmml::INF, inputOptions.gapCTermination_ + " cap for : " + residue->printId());
                // Do it
                chain->InsertCap(*previous, inputOptions.gapCTermination_);
                chain->InsertCap(*residue, inputOptions.gapNTermination_);
                // Record it
                ppInfo.missingResidues_.emplace_back(previous->getChainId(), previous->getNumberAndInsertionCode(), residue->getNumberAndInsertionCode(), inputOptions.gapCTermination_, inputOptions.gapNTermination_);
            }
            previous = residue;
            previousSequenceNumber = residue->getNumber();
        }
    }
}

void PdbModel::preProcessMissingUnrecognized(pdb::PreprocessorInformation &ppInfo)
{
    parameters::Manager parmManager;
    for(auto &residue : this->getResidues())
    {
        std::vector<std::string> parmAtomNames = parmManager.GetAtomNamesForResidue(residue->GetParmName());
        std::vector<std::string> parmHeavyAtomNames = parmManager.GetHeavyAtomNamesForResidue(residue->GetParmName());
        // Unrecognized residue->
        if (parmAtomNames.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "ParmManager did not recognize residue: " + residue->GetParmName());
            ppInfo.unrecognizedResidues_.emplace_back(residue->getId());
        }
        else // Recognized residue->
        {
            std::vector<std::string> pdbAtomNames = residue->getAtomNames();
            for (auto &parmHeavyAtomName : parmHeavyAtomNames) // What heavy atoms are missing from the pdb residue?
            {
                if ( std::find(pdbAtomNames.begin(), pdbAtomNames.end(), parmHeavyAtomName) == pdbAtomNames.end() )
                { // Residue missing a heavy atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Atom named " + parmHeavyAtomName + " missing from " + residue->printId());
                    ppInfo.missingHeavyAtoms_.emplace_back(parmHeavyAtomName, residue->getId());
                }
            }
            for (auto &pdbAtomName : pdbAtomNames) // What atoms in the pdb residue are unrecognized?
            {
                if ( std::find(parmAtomNames.begin(), parmAtomNames.end(), pdbAtomName) == parmAtomNames.end() )
                {
                    // Residue contains unrecognized atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Unrecognized atom named " + pdbAtomName + " in " + residue->printId());
                    ppInfo.unrecognizedAtoms_.emplace_back(pdbAtomName, residue->getId());
                }
            }
        }
    }
    return;
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////
void PdbModel::Print(std::ostream &out) const
{
    for (auto &residue : this->getResidues())
    {
        residue->Print(out);
    }
}

void PdbModel::Write(std::ostream& stream) const
{
    stream << "MODEL " << std::right << std::setw(4) << this->getNumber() << "\n";
    for (auto &chain : this->getMolecules())
    {
        chain->Write(stream);
    }
    for (auto &conect : this->GetConectRecords())
    {
        conect.Write(stream);
    }
    stream << "ENDMDL\n";
    return;
}

