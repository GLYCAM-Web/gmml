#include <algorithm> // std::find
#include "includes/InputSet/PdbFile/coordinateSection.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/common.hpp" // gmml::PROTEINS
using pdb::CoordinateSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CoordinateSection::CoordinateSection() {}

CoordinateSection::CoordinateSection(std::stringstream &stream_block)
{
    int currentModelNumber = 1;
//    PdbResidue* currentResidue;
    std::vector<std::unique_ptr<AtomRecord>> atomRecords;
    std::string currentResidueId = "?_?_?_?_?"; // model_resname_insertionCode_resnumber_chain_
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
            atomRecords_.push_back(std::make_unique<AtomRecord>(line, currentModelNumber));
            atomRecords.push_back(std::make_unique<AtomRecord>(line, currentModelNumber));
            if ( (previousResidueId != "InitialValue") && (previousResidueId != atomRecords.back()->GetResidueId()) )
            {
                //residues_.push_back(std::make_unique<PdbResidue>(atomRecords));
                residues_.emplace_back(atomRecords);
            }
//            residues_.emplace_back(currentModelNumber, line, stream_block);

 //           atomRecords_.emplace_back(line, currentModelNumber);
            // This next thing is to pre-organize into residues. Maybe dumb.
//            currentResidueId = newAtomRecord.GetResidueId();
//            if(currentResidueId == previousResidueId)
//            {
//                currentResidue->AddAtom(&newAtomRecord);
//                std::cout << "Added to current residue\n";
//            }
//            else
//            {
//                currentResidue = this->CreateNewResidue(&newAtomRecord);
//                std::cout << "New residue created!\n";
//            }
//            previousResidueId = currentResidueId;
        }
        std::cerr << "Residues so far are:\n";
        for (auto &residue : residues_)
        {
            residue->Print();
        }
    }
}


std::vector<std::vector<pdb::PdbResidue>> CoordinateSection::GetModels() const
{
    std::vector<std::vector<pdb::PdbResidue>> models;
    int previousModel = -123456789; // It would need to be very wrong to be this value.
    for(auto &residue : this->GetResidues())
    {
        if (residue.GetModelNumber() != previousModel )
        {
            models.emplace_back(std::vector<PdbResidue>{residue});
        }
        else
        {
            models.back().push_back(residue);
        }
        previousModel = residue.GetModelNumber();
    }
    return models;
}

std::vector<std::vector<pdb::PdbResidue>> CoordinateSection::GetProteinChains()
{
    std::vector<std::vector<pdb::PdbResidue>> chains;
    std::string previousChain = "AUniqueLittleString";
    int previousModelNumber = -99999;
    for(auto residue : this->GetResidues())
    {
        // This gmml::PROTEINS seems weird, but whatever, it works.
        if( std::find( gmml::PROTEINS, ( gmml::PROTEINS + gmml::PROTEINSSIZE ), residue.GetName() ) != ( gmml::PROTEINS + gmml::PROTEINSSIZE ) )
        {
            if ( residue.GetChainId() != previousChain || residue.GetModelNumber() != previousModelNumber )
            {
                chains.emplace_back(std::vector<PdbResidue>{residue});
            }
            else
            {
                chains.back().push_back(residue);
            }
        }
        previousChain = residue.GetChainId();
        previousModelNumber = residue.GetModelNumber();
    }
    // Labels
    for(auto &chain : chains)
    {
        chain.front().AddLabel("NTerminal");
        chain.back().AddLabel("CTerminal");
    }
    return chains;
}

std::vector<pdb::PdbResidue> CoordinateSection::GetResidues() const
{
    std::vector<pdb::PdbResidue> residues;
    std::string id = "?_?_?_?_?"; // resname_insertionCode_resnumber_chain_model
    std::string previousId = "AStringThatIsUniqueFromIdForTheFirstLoop";
    for(auto &atomUniquePtr : atomRecords_)
    {
        id = atomUniquePtr->GetResidueId();
        if(id == previousId)
        {
            residues.back().AddAtom(atomUniquePtr.get());
        }
        else
        {
            residues.emplace_back(atomUniquePtr.get());
        }
        previousId = id;
    }
    return residues;
}

std::vector<pdb::PdbResidue> CoordinateSection::FindResidues(const std::string selector)
{
    std::vector<pdb::PdbResidue> matchingResidues;
    for(auto &residue : this->GetResidues())
    {
        std::size_t found = residue.GetId().find(selector);
        if(found != std::string::npos)
        {
            matchingResidues.push_back(residue);
            //std::cout << residue.GetId() << " contains " << selector << "\n";
        }
    }
    return matchingResidues;
}

void CoordinateSection::ChangeResidueName(const std::string& selector, const std::string& newName)
{
    for(auto &residue : this->GetResidues())
    {
        std::size_t found = residue.GetId().find(selector);
        if(found != std::string::npos)
        {
            residue.SetName(newName);
        }
    }
    return;
}

pdb::AtomRecord* CoordinateSection::FindAtom(int serialNumber)
{
    AtomRecord* nullRecord = nullptr;
    for(auto &atomRecord : atomRecords_)
    {
        if (atomRecord->GetSerialNumber() == serialNumber )
        {
            return atomRecord.get(); // get() gets the raw ptr, imo raw() would have been cooler and clearer.
        }
    }
    return nullRecord;
}

void CoordinateSection::DeleteAtomRecord(AtomRecord* atom)
{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
    auto i = this->FindPositionOfAtom(atom);
    if (i != atomRecords_.end())
    {
       i = atomRecords_.erase(i);
       gmml::log(__LINE__,__FILE__,gmml::INF, "Atom " + atom->GetId() + " has been erased. You're welcome.");
    }
    else
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not delete " + atom->GetId() + " as it was not found in atom records\n");
    }
    return;
}

pdb::AtomRecordIterator CoordinateSection::CreateNewAtomRecord(std::string name, GeometryTopology::Coordinate& coord, AtomRecord* sisterAtom)
{ // A bit wonky as position records the sisterAtom position in the vector, and then the new atom position to be returned.
    AtomRecordIterator position = this->FindPositionOfAtom(sisterAtom);
    if (position != atomRecords_.end())
    {
        ++position; // it is ok to insert at end(). I checked. It was ok. Ok.
        position = atomRecords_.insert(position, std::make_unique<AtomRecord>(name, coord, sisterAtom));
        gmml::log(__LINE__,__FILE__,gmml::INF, "New atom named " + name + " has been born; You're welcome.");
    }
    else
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create atom named " + name + " as sisterAtom was not found in atom records\n");
    }
    return position;
}

pdb::AtomRecordIterator CoordinateSection::CreateNewAtomRecord(const std::string& atomName, const std::string& residueName, const int& residueSequenceNumber, const GeometryTopology::Coordinate& coord, const std::string& chainId, const int& modelNumber, AtomRecordIterator previousAtomPosition)
{
    //auto position = this->FindPositionOfAtom(previousAtom);
    //if (position != atomRecords_.end())
    //{
        AtomRecordIterator newAtomPosition = atomRecords_.insert(++previousAtomPosition, std::make_unique<AtomRecord>(atomName, residueName, residueSequenceNumber, coord, chainId, modelNumber));
        gmml::log(__LINE__,__FILE__,gmml::INF, "New atom named " + atomName + " has been born; You're welcome.");
    //}
//    else
//    {
//        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create atom named " + atomName + " as previousAtom was not found in atom records\n");
//    }
    return newAtomPosition; // this is an iterator to where the unique ptr is, and get() returns a raw ptr.
}


std::vector<std::unique_ptr<pdb::AtomRecord>>::iterator CoordinateSection::FindPositionOfAtom(AtomRecord* queryAtom)
{
    auto i = atomRecords_.begin();
    auto e = atomRecords_.end();
    while (i != e)
    {
        if (queryAtom == i->get())
        {
            return i;
        }
        else
        {
            ++i;
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find " + queryAtom->GetId() + " in atom records\n");
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
void CoordinateSection::Print(std::ostream &out) const
{
    out << "The atom records are: " << "\n";
    for (auto &atomRecord : atomRecords_)
    {
        atomRecord->Print(out);
    }
}
void CoordinateSection::Write(std::ostream& stream) const
{
    std::vector<std::vector<pdb::PdbResidue>> models = this->GetModels(); // Vectors of residues in a vector organized by model.
    for (auto &model : models)
    {
        if (models.size() > 1)
        {
            stream << "MODEL " << std::right << std::setw(4) << model.at(0).GetModelNumber() << "\n";
        }
        std::string previousResidueName = "";
        std::string previousChainId = model.at(0).GetChainId();
        for (auto &residue: model) // Maybe I should just put the TER cards in the vector as AtomRecord?
        { // if previous was NME, or it's a new chain, or it's a HETATM and starts previous residue name began a 0 (e.g. glycam 0SA), insert a TER.
            if (previousResidueName == "NME" || previousChainId != residue.GetChainId() || (residue.GetRecordName() == "HETATM" && previousResidueName.substr(0,1) == "0"))
            {
                stream << "TER\n";
            }
            residue.Write(stream);
            previousResidueName = residue.GetName();
            previousChainId = residue.GetChainId();
        }
        if (models.size() > 1)
        {
            stream << "ENDMDL\n";
        }
    }
}

