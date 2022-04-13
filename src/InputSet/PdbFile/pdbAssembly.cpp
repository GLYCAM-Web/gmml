#include <includes/InputSet/PdbFile/pdbAssembly.hpp>
#include <algorithm> // std::find
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/common.hpp" // gmml::PROTEINS
using pdb::PdbAssembly;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAssembly::PdbAssembly() {}

PdbAssembly::PdbAssembly(std::stringstream &stream_block)
{
    int currentModelNumber = 1;
//    PdbResidue* currentResidue;
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
                gmml::log(__LINE__, __FILE__, gmml::ERR, "Model issue: this ain't an int: " + codeUtils::RemoveWhiteSpace(line.substr(10,4)));
                currentModelNumber = 1; // Seems like a reasonable default.
            }
        }
        // ATOM
        else if ( (recordName == "ATOM") || (recordName == "HETATM") )
        {
//            atomRecords_.push_back(std::make_unique<AtomRecord>(line, currentModelNumber));
            std::string residueId = this->PeekAtResidueId(line);
            if (previousResidueId != residueId)
            { // Create a residue;
                residues_.push_back(std::make_unique<PdbResidue>(line, currentModelNumber));
                previousResidueId = residueId;
            }
            else
            { // Add to previous residue;
                residues_.back()->CreateAtom(line, currentModelNumber);
            }
        }
        // TER . Indicating break in chain i.e. start of new molecule or a branch
        else if ( (recordName == "TER") && (!residues_.empty()) )
        {
            residues_.back()->AddTerCard();
        }
    }
//    for (auto &residue : residues_)
//    {
//        residue->Print();
//    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

// This is bad as it repeats how to read a line and how to create a residue ID, but I need to know which residue to put the atom into before I construct the atom.
std::string PdbAssembly::PeekAtResidueId(const std::string &line)
{
    // Dealing with number overruns for serialNumber and residueNumber
    int shift = codeUtils::GetSizeOfIntInString(line.substr(12));
    std::string residueName = codeUtils::RemoveWhiteSpace(line.substr(17 + shift, 3));
    std::string chainId = codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
    int secondShift = codeUtils::GetSizeOfIntInString(line.substr(26 + shift));
    std::string residueNumber = codeUtils::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift));
    // Insertion code gets shifted right by every overrun in residue number.
    std::string insertionCode = codeUtils::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
    return residueName + "_" + residueNumber + "_" + insertionCode + "_" + chainId;
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

std::vector<std::vector<pdb::PdbResidue*>> PdbAssembly::GetProteinChains()
{
    std::vector<std::vector<pdb::PdbResidue*>> chains;
    std::string previousChain = "AUniqueLittleString";
    int previousModelNumber = -99999;
    for(auto &residue : residues_)
    {
        // This gmml::PROTEINS seems weird, but whatever, it works.
        if( std::find( gmml::PROTEINS, ( gmml::PROTEINS + gmml::PROTEINSSIZE ), residue->GetName() ) != ( gmml::PROTEINS + gmml::PROTEINSSIZE ) )
        {
            if ( residue->GetChainId() != previousChain || residue->GetModelNumber() != previousModelNumber )
            {
                chains.emplace_back(std::vector<PdbResidue*>{residue.get()});
                previousChain = residue->GetChainId();
                previousModelNumber = residue->GetModelNumber();
            }
            else
            {
                chains.back().push_back(residue.get());
            }
        }
    }
    // Labels
    for(auto &chain : chains)
    {
        chain.front()->addLabel("NTerminal");
        chain.back()->addLabel("CTerminal");
    }
    return chains;
}

//std::vector<pdb::PdbResidue*> PdbAssembly::GetResidues() const
//{
//    std::vector<pdb::PdbResidue*> residues; // reserve the right amount of memory.
//    for(auto &resUniquePtr : residues_)
//    {
//        residues.push_back(resUniquePtr.get());
//    }
//    return residues;
//}

std::vector<pdb::PdbResidue*> PdbAssembly::FindResidues(const std::string selector)
{
    std::vector<pdb::PdbResidue*> matchingResidues;
    for(auto &residue : residues_)
    {
        std::size_t found = residue->GetId().find(selector);
        if(found != std::string::npos)
        {
            matchingResidues.push_back(residue.get());
        }
    }
    return matchingResidues;
}

void PdbAssembly::ChangeResidueName(const std::string& selector, const std::string& newName)
{
    for(auto &residue : residues_)
    {
        std::size_t found = residue->GetId().find(selector);
        if(found != std::string::npos)
        {
            residue->SetName(newName);
            return;
        }
    }
    return;
}

pdb::AtomRecord* PdbAssembly::FindAtom(const int& serialNumber) const
{
    AtomRecord* foundAtom = nullptr;
    for(auto &residue : residues_)
    {
        foundAtom = residue->FindAtom(serialNumber);
        if ( foundAtom != nullptr )
        {
            return foundAtom; // get() gets the raw ptr, imo raw() would have been cooler and clearer.
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


pdb::PdbResidue* PdbAssembly::CreateNewResidue(const std::string residueName, const std::string atomName, GeometryTopology::Coordinate& atomCoord, const pdb::PdbResidue& referenceResidue)
{
    //Where the residue is in the vector matters. It should go after the reference residue.
    pdb::PdbResidueIterator position = this->FindPositionOfResidue(&referenceResidue);
    if (position != residues_.end())
    {
        ++position; // it is ok to insert at end(). I checked. It was ok. Ok.
        position = residues_.insert(position, std::make_unique<PdbResidue>(residueName, atomName, atomCoord, &referenceResidue));
        gmml::log(__LINE__,__FILE__,gmml::INF, "New residue named " + residueName + " has been born; You're welcome.");
    }
    else
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create residue named " + residueName + " as referenceResidue was not found\n");
    }
    return (*position).get(); // Wow ok, so dereference the reference to a uniquePtr, then use get() to create a raw ptr.
}

pdb::PdbResidueIterator PdbAssembly::FindPositionOfResidue(const PdbResidue* queryResidue)
{
    auto i = residues_.begin();
    auto e = residues_.end();
    while (i != e)
    {
        if (queryResidue == i->get())
        {
            return i;
        }
        else
        {
            ++i;
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find " + queryResidue->GetId() + " in atom records\n");
    return e;
}
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
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////
void PdbAssembly::Print(std::ostream &out) const
{
    for (auto &residue : residues_)
    {
        residue->Print(out);
    }
}
void PdbAssembly::Write(std::ostream& stream) const
{
    std::vector<std::vector<pdb::PdbResidue*>> models = this->GetModels(); // Vectors of residues in a vector organized by model.
    for (auto &model : models)
    {
        if (models.size() > 1)
        {
            stream << "MODEL " << std::right << std::setw(4) << model.at(0)->GetModelNumber() << "\n";
        }
        std::string previousResidueName = "";
        std::string previousChainId = model.at(0)->GetChainId();
        for (auto &residue: model) // Maybe I should just put the TER cards in the vector as AtomRecord?
        { // if previous was NME, or it's a new chain, or it's a HETATM and starts previous residue name began a 0 (e.g. glycam 0SA), insert a TER.
            if (previousChainId != residue->GetChainId() || (residue->GetRecordName() == "HETATM" && previousResidueName.substr(0,1) == "0"))
            {
                stream << "TER\n";
            }
            residue->Write(stream);
            previousResidueName = residue->GetName();
            previousChainId = residue->GetChainId();
        }
        if (models.size() > 1)
        {
            stream << "ENDMDL\n";
        }
    }
}

