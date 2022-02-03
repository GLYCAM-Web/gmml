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
//    std::string currentResidueId = "?_?_?_?_?"; // model_resname_insertionCode_resnumber_chain_
//    std::string previousResidueId = "IJustNeedThisToBeDifferentFromCurrentResidueIdForTheinitialCheckOkDontJudgeMeImDoingMyBestOutHere";
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
 //           atomRecords_.emplace_back(line, currentModelNumber);
            atomRecords_.push_back(std::make_unique<AtomRecord>(line, currentModelNumber));
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
//        std::cerr << "Residues so far are:\n";
//        for (auto &residue : residues_)
//        {
//            residue.Print();
//        }
    }
}
std::vector<std::vector<pdb::PdbResidue>> CoordinateSection::GetProteinChains()
{
    std::vector<std::vector<pdb::PdbResidue>> chains;
    std::string previousChain = "AUniqueLittleString";
    for(auto residue : this->GetResidues())
    {
        // This gmml::PROTEINS seems weird, but whatever, it works.
        if( std::find( gmml::PROTEINS, ( gmml::PROTEINS + gmml::PROTEINSSIZE ), residue.GetName() ) != ( gmml::PROTEINS + gmml::PROTEINSSIZE ) )
        {
            if (residue.GetChainId() == previousChain)
            {
                chains.back().push_back(residue);
            }
            else
            {
                chains.emplace_back(std::vector<PdbResidue>{residue});
            }
        }
        previousChain = residue.GetChainId();
    }
    // Labels
    for(auto &chain : chains)
    {
        chain.front().addLabel("NTerminal");
        chain.back().addLabel("CTerminal");
    }
    return chains;
}

std::vector<pdb::PdbResidue> CoordinateSection::GetResidues()
{
    std::vector<pdb::PdbResidue> residues;
    std::string id = "?_?_?_?_?"; // model_resname_insertionCode_resnumber_chain_
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
            //std::cout << "Creating new residue\n";
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

void CoordinateSection::CreateNewAtomRecord(std::string name, GeometryTopology::Coordinate& coord, AtomRecord* sisterAtom)
{
    auto position = this->FindPositionOfAtom(sisterAtom);
    if (position != atomRecords_.end())
    {
        ++position; // it is ok to insert at end(). I checked. It was ok. Ok.
        atomRecords_.insert(position, std::make_unique<AtomRecord>(name, coord, sisterAtom));
        gmml::log(__LINE__,__FILE__,gmml::INF, "New atom named " + name + " has been born; You're welcome.");
    }
    else
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create atom named " + name + " as sisterAtom was not found in atom records\n");
    }
    return;
}

pdb::AtomRecord* CoordinateSection::CreateNewAtomRecord(const std::string& atomName, const std::string& residueName, const int& residueSequenceNumber, const GeometryTopology::Coordinate& coord, const std::string& chainId, const int& modelNumber, AtomRecord* previousAtom)
{
    auto position = this->FindPositionOfAtom(previousAtom);
    if (position != atomRecords_.end())
    {
        ++position;
        position = atomRecords_.insert(position, std::make_unique<AtomRecord>(atomName, residueName, residueSequenceNumber, coord, chainId, modelNumber));
        gmml::log(__LINE__,__FILE__,gmml::INF, "New atom named " + atomName + " has been born; You're welcome.");
    }
    else
    {
        gmml::log(__LINE__,__FILE__,gmml::ERR, "Could not create atom named " + atomName + " as previousAtom was not found in atom records\n");
    }
    return position->get(); // this is an iterator to where the unique ptr is, and get() returns a raw ptr.
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
//    out << "The number of residues is: " << this->GetResidues().size() << "\n";
//    out << "The residues are: " << "\n";
//    for (auto &residue : this->GetResidues())
//    {
//        residue.Print();
//    }
}
