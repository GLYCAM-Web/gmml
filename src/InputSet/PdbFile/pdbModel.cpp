#include <algorithm> // std::find

#include "includes/InputSet/PdbFile/pdbModel.hpp"
#include "includes/InputSet/PdbFile/pdbChain.hpp"
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/common.hpp" // gmml::PROTEINS

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
        }
        // ATOM
        else if ( (recordName == "ATOM") || (recordName == "HETATM") )
        {
            // Gimme everything with the same chain, can be everything with no chain.
            // Function that will read from stringstream until chain ID changes or TER or just not ATOM/HETATM
            std::stringstream singleChainSection = this->extractSingleChainFromRecordSection(stream_block, line, this->extractChainId(line));
            //this->addMolecule(PdbChain(singleChainSection, this->extractChainId(line)));
            gmml::log(__LINE__,__FILE__,gmml::INF, "I'm about to die now?");
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

std::stringstream PdbModel::extractSingleChainFromRecordSection(std::stringstream &pdbFileStream, std::string line, const std::string& initialChainID)
{
    std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
    std::stringstream singleChainSection;
    std::string chainID = initialChainID;
    while(chainID == initialChainID)
    {
        singleChainSection << line << std::endl;
        previousLinePosition = pdbFileStream.tellg(); // Save current line position.
        if(!std::getline(pdbFileStream, line))
        {
            break; // // If we hit the end, time to leave.
        }
        chainID = this->extractChainId(line);
    }
    pdbFileStream.seekg(previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    gmml::log(__LINE__,__FILE__,gmml::INF, singleChainSection.str());
    return singleChainSection;
}

void PdbModel::addConectRecord(const AtomRecord* atom1, const AtomRecord* atom2)
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

//std::vector<std::vector<pdb::PdbResidue*>> PdbModel::GetProteinChains()
//{
//    std::vector<std::vector<pdb::PdbResidue*>> chains;
//    std::string previousChain = "AUniqueLittleString";
//    int previousModelNumber = -99999;
//    for(auto &residue : this->getResidues())
//    {
//        // This gmml::PROTEINS seems weird, but whatever, it works.
//        if( std::find( gmml::PROTEINS, ( gmml::PROTEINS + gmml::PROTEINSSIZE ), residue->GetName() ) != ( gmml::PROTEINS + gmml::PROTEINSSIZE ) )
//        {
//            if ( residue->GetChainId() != previousChain || residue->GetModelNumber() != previousModelNumber )
//            {
//                chains.emplace_back(std::vector<PdbResidue*>{residue.get()});
//                previousChain = residue->GetChainId();
//                previousModelNumber = residue->GetModelNumber();
//            }
//            else
//            {
//                chains.back().push_back(residue.get());
//            }
//        }
//    }
//    // Labels
//    for(auto &chain : chains)
//    {
//        chain.front()->addLabel("NTerminal");
//        chain.back()->addLabel("CTerminal");
//    }
//    return chains;
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
//void PdbModel::ChangeResidueName(const std::string& selector, const std::string& newName)
//{
//    for(auto &residue : residues_)
//    {
//        std::size_t found = residue->GetId().find(selector);
//        if(found != std::string::npos)
//        {
//            residue->SetName(newName);
//            return;
//        }
//    }
//    return;
//}
//
const pdb::AtomRecord* PdbModel::FindAtom(const int& serialNumber) const
{
    for(auto &atom : this->getAtoms())
    {
        if (atom->GetSerialNumber() == serialNumber)
        {
            return atom;
        }
    }
    return nullptr;
}

//void CoordinateSection::DeleteAtomRecord(AtomRecord* atom)
//{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
//    for(auto &residue : residues_)
//    {
//        if (residue->DeleteAtomRecord(atom)) // returns true if it found and deleted the atomRecord
//        {
//            return;
//        }
//    }
//    return;
//}

//pdb::AtomRecordIterator CoordinateSection::CreateNewAtomRecord(std::string name, GeometryTopology::Coordinate& coord, AtomRecord* sisterAtom)
//{ // A bit wonky as position records the sisterAtom position in the vector, and then the new atom position to be returned.
//    AtomRecordIterator position = this->FindPositionOfAtom(sisterAtom);
//    if (position != atomRecords_.end())
//    {
//        ++position; // it is ok to insert at end(). I checked. It was ok. Ok.
//        position = atomRecords_.insert(position, std::make_unique<AtomRecord>(name, coord, sisterAtom));
//        gmml::log(__LINE__,__FILE__,gmml::INF, "New atom named " + name + " has been born; You're welcome.");
//    }
//    else
//    {
//        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create atom named " + name + " as sisterAtom was not found in atom records\n");
//    }
//    return position;
//}

//pdb::AtomRecordIterator CoordinateSection::CreateNewAtomRecord(const std::string& atomName, const std::string& residueName, const int& residueSequenceNumber, const GeometryTopology::Coordinate& coord, const std::string& chainId, const int& modelNumber, AtomRecordIterator previousAtomPosition)
//{
//    //auto position = this->FindPositionOfAtom(previousAtom);
//    //if (position != atomRecords_.end())
//    //{
//        AtomRecordIterator newAtomPosition = atomRecords_.insert(++previousAtomPosition, std::make_unique<AtomRecord>(atomName, residueName, residueSequenceNumber, coord, chainId, modelNumber));
//        gmml::log(__LINE__,__FILE__,gmml::INF, "New atom named " + atomName + " has been born; You're welcome.");
//    //}
////    else
////    {
////        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create atom named " + atomName + " as previousAtom was not found in atom records\n");
////    }
//    return newAtomPosition; // this is an iterator to where the unique ptr is, and get() returns a raw ptr.
//}


//pdb::PdbResidue* PdbModel::CreateNewResidue(const std::string residueName, const std::string atomName, GeometryTopology::Coordinate& atomCoord, const pdb::PdbResidue& referenceResidue)
//{
//    //Where the residue is in the vector matters. It should go after the reference residue.
//    pdb::PdbResidueIterator position = this->FindPositionOfResidue(&referenceResidue);
//    if (position != residues_.end())
//    {
//        ++position; // it is ok to insert at end(). I checked. It was ok. Ok.
//        position = residues_.insert(position, std::make_unique<PdbResidue>(residueName, atomName, atomCoord, &referenceResidue));
//        gmml::log(__LINE__,__FILE__,gmml::INF, "New residue named " + residueName + " has been born; You're welcome.");
//    }
//    else
//    {
//        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create residue named " + residueName + " as referenceResidue was not found\n");
//    }
//    return (*position).get(); // Wow ok, so dereference the reference to a uniquePtr, then use get() to create a raw ptr.
//}
//
//
//pdb::PdbResidueIterator PdbModel::FindPositionOfResidue(const PdbResidue* queryResidue)
//{
//    auto i = residues_.begin();
//    auto e = residues_.end();
//    while (i != e)
//    {
//        if (queryResidue == i->get())
//        {
//            return i;
//        }
//        else
//        {
//            ++i;
//        }
//    }
//    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find " + queryResidue->GetId() + " in atom records\n");
//    return e;
//}
//////////////////////////////////////////////////////////
//                      ACCESSOR                        //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      MUTATOR                         //
//////////////////////////////////////////////////////////
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
        AtomRecord* sgAtom1 = cysRes1->FindAtom("SG");
        for (std::vector<pdb::PdbResidue*>::iterator it2 = std::next(it1, 1); it2 != cysResidues.end(); ++it2)
        {
            PdbResidue* cysRes2 = *it2;
            AtomRecord* sgAtom2 = cysRes2->FindAtom("SG");
            if ( (sgAtom1 != nullptr) && (sgAtom2 != nullptr) )
            {
                //gmml::log(__LINE__, __FILE__, gmml::INF, "Found SG ATOMS");
                double distance = sgAtom1->CalculateDistance(sgAtom2);
                if (distance < gmml::dSulfurCutoff && distance > 0.001)
                {
                  //  gmml::log(__LINE__, __FILE__, gmml::INF, "Distance less than cutoff");
                    cysRes1->setName("CYX");
                    cysRes2->setName("CYX");
                   // gmml::log(__LINE__, __FILE__, gmml::INF, "Names set");
                    this->addConectRecord(sgAtom1, sgAtom2);
                    ppInfo.cysBondResidues_.emplace_back(cysRes1, cysRes2, distance);
                   // gmml::log(__LINE__, __FILE__, gmml::INF, "ThisNoHappen?");
                    std::stringstream message;
                    message << "Bonding " << cysRes1->getId() << " and " << cysRes2->getId() << " with distance " << distance;
                    gmml::log(__LINE__, __FILE__, gmml::INF, message.str());
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////
void PdbModel::Print(std::ostream &out) const
{
    for (auto &residue : this->getResidues())
    {
       // residue->Print(out);
    }
}

void PdbModel::Write(std::ostream& stream) const
{
    stream << "MODEL " << std::right << std::setw(4) << this->getModelNumber() << "\n";
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

