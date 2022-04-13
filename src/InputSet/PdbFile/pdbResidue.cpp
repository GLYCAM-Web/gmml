#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/CodeUtils/logging.hpp"
using pdb::PdbResidue;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
//PdbResidue::PdbResidue(AtomRecord* atomRecord)
//{
//    this->AddAtom(atomRecord);
//    hasTerCard_ = false;
//}
//
//PdbResidue::PdbResidue(std::vector<AtomRecord*> atomRecords)
//{
//    atoms_ = atomRecords;
//    hasTerCard_ = false;
//}
PdbResidue::PdbResidue(const std::string& line, const int& modelNumber)
{
    //atoms_.push_back(std::make_unique<AtomRecord>(line, modelNumber));
    this->CreateAtomFromLine(line, modelNumber);
    modelNumber_ = modelNumber;
    hasTerCard_ = false;
}

PdbResidue::PdbResidue(const std::string residueName, const std::string atomName, GeometryTopology::Coordinate& atomCoord, const PdbResidue *referenceResidue)
{
    AtomRecord tempAtom(atomName, residueName, referenceResidue->GetSequenceNumber() + 1, referenceResidue->GetInsertionCode(), atomCoord, referenceResidue->GetChainId(), referenceResidue->GetModelNumber());
    this->createAtom(tempAtom); // This should call the cdsResidue CreateAtom function.
    //atoms_.push_back(std::make_unique<AtomRecord>(atomName, residueName, referenceResidue->GetSequenceNumber() + 1, referenceResidue->GetInsertionCode(), atomCoord, referenceResidue->GetChainId(), referenceResidue->GetModelNumber()));
    modelNumber_ = referenceResidue->GetModelNumber();
    hasTerCard_ = false;
}

//PdbResidue::PdbResidue(std::vector<std::unique_ptr<AtomRecord>>& atomRecords)
//{
//    atomRecordss_.swap(atomRecords); // atomRecords will become empty, atoms_ will contain contents of atomRecords.
//    hasTerCard_ = false;
//}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
//std::vector<std::string> PdbResidue::GetAtomNames() const
//{
//    std::vector<std::string> foundAtomNames;
//    for(auto &atomRecord : this->getAtoms())
//    {
//        foundAtomNames.push_back(atomRecord->GetName());
//    }
//    return foundAtomNames;
//}

//pdb::AtomRecord* PdbResidue::GetLastAtom() const
//{
//    if (atoms_.empty())
//    {
//        gmml::log(__LINE__, __FILE__, gmml::ERR, "atomRecords in PdbResidue is empty. Impossible!?");
//    }
//    return atoms_.back().get();
//}
//
//pdb::AtomRecord* PdbResidue::GetFirstAtom() const
//{
//    if (atoms_.empty())
//    {
//        gmml::log(__LINE__, __FILE__, gmml::ERR, "atomRecords in PdbResidue is empty. Impossible!?");
//    }
//    return atoms_.front().get();
//}

//const std::string& PdbResidue::GetChainId() const
//{
//    return this->GetFirstAtom()->GetChainId();
//}
//const std::string& PdbResidue::GetName() const
//{
//    return this->GetFirstAtom()->GetResidueName();
//}

//const std::string& PdbResidue::GetRecordName() const
//{
//    return this->GetFirstAtom()->GetRecordName();
//}

const std::string PdbResidue::GetParmName() const // If terminal, need to look up e.g. NPRO or CPRO instead of PRO.
{
    if (this->containsLabel("NTerminal"))
    {
        return "N" + this->GetName();
    }
    else if (this->containsLabel("CTerminal"))
    {
        return "C" + this->GetName();
    }
    return this->GetName();
}
//const int& PdbResidue::GetSequenceNumber() const
//{
//    return this->GetFirstAtom()->GetResidueSequenceNumber();
//}

//const std::string& PdbResidue::GetInsertionCode() const
//{
//    return this->GetFirstAtom()->GetInsertionCode();
//}
//
//std::string PdbResidue::GetId() const
//{
//    return this->GetFirstAtom()->GetResidueId();
//}

//const int& PdbResidue::GetModelNumber() const
//{
//    return modelNumber_;
//}

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidue::CreateAtomFromLine(const std::string& line, const int& currentModelNumber)
{
    AtomRecord tempRecord(line, currentModelNumber);
    this->createAtom(tempRecord);
    //atoms_.push_back(std::make_unique<AtomRecord>(line, currentModelNumber));
    return;
}

//void PdbResidue::CreateAtom(const std::string atomName, GeometryTopology::Coordinate& atomCoord)
//{
//    atoms_.push_back(std::make_unique<AtomRecord>(atomName, this->GetName(), this->GetSequenceNumber(), this->GetInsertionCode(), atomCoord, this->GetChainId(), this->GetModelNumber()));
//    return;
//}


//void PdbResidue::AddAtom(AtomRecord* atomRecord)
//{
//    //std::cout << "Adding atom with id " << atomRecord->GetId() << " to residue.\n";
//    atoms_.push_back(atomRecord);
//    return;
//}

//void PdbResidue::SetName(const std::string name)
//{
//    for(auto &atom : atoms_)
//    {
//        atom->SetResidueName(name);
//    }
//    return;
//}

//bool PdbResidue::DeleteAtomRecord(AtomRecord* atom)
//{ // Passing in a raw ptr, but the vector is unique_ptr so gotta use i->get() to compare raws.
//    auto i = this->FindPositionOfAtom(atom);
//    if (i != atoms_.end())
//    {
//       i = atoms_.erase(i);
//       gmml::log(__LINE__,__FILE__,gmml::INF, "Atom " + atom->GetId() + " has been erased. You're welcome.");
//       return true;
//    }
//    return false;
//}

//std::vector<std::unique_ptr<pdb::AtomRecord>>::iterator PdbResidue::FindPositionOfAtom(AtomRecord* queryAtom)
//{
//    auto i = atoms_.begin();
//    auto e = atoms_.end();
//    while (i != e)
//    {
//        if (queryAtom == i->get())
//        {
//            return i;
//        }
//        else
//        {
//            ++i;
//        }
//    }
//    gmml::log(__LINE__,__FILE__,gmml::ERR, "Did not find " + queryAtom->GetId() + " in atom records\n");
//    return e;
//}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidue::Print(std::ostream &out) const
{
    out << "pdb::Residue : " << this->GetId() << std::endl;
    for(auto &atom : this->getAtoms())
    {
        out << "    atom : " << atom->GetId() << " X: "  << atom->GetCoordinate().GetX() << " Y: "  << atom->GetCoordinate().GetY() << " Z: " << atom->GetCoordinate().GetZ() << "\n";
    }
}

void PdbResidue::Write(std::ostream& stream) const
{
    for(auto &atom : this->getAtoms())
    {
        atom->Write(stream);
    }
    if(hasTerCard_)
    {
        stream << "TER\n";
    }
}
