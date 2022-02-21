#include "includes/InputSet/PdbFile/pdbResidue.hpp"
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
//    atomRecords_ = atomRecords;
//    hasTerCard_ = false;
//}
PdbResidue::PdbResidue(const std::string& line, const int& modelNumber)
{
    atomRecords_.push_back(std::make_unique<AtomRecord>(line, modelNumber));
    modelNumber_ = modelNumber;
    hasTerCard_ = false;
}
//PdbResidue::PdbResidue(std::vector<std::unique_ptr<AtomRecord>>& atomRecords)
//{
//    atomRecordss_.swap(atomRecords); // atomRecords will become empty, atomRecords_ will contain contents of atomRecords.
//    hasTerCard_ = false;
//}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::vector<std::string> PdbResidue::GetAtomNames() const
{
    std::vector<std::string> foundAtomNames;
    for(auto &atomRecord : atomRecords_)
    {
        foundAtomNames.push_back(atomRecord->GetName());
    }
    return foundAtomNames;
}

pdb::AtomRecord* PdbResidue::GetLastAtom() const
{
    if (atomRecords_.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "atomRecords in PdbResidue is empty. Impossible!?");
    }
    return atomRecords_.back().get();
}

pdb::AtomRecord* PdbResidue::GetFirstAtom() const
{
    if (atomRecords_.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "atomRecords in PdbResidue is empty. Impossible!?");
    }
    return atomRecords_.front().get();
}

const std::string& PdbResidue::GetChainId() const
{
    return this->GetFirstAtom()->GetChainId();
}
const std::string& PdbResidue::GetName() const
{
    return this->GetFirstAtom()->GetResidueName();
}

const std::string& PdbResidue::GetRecordName() const
{
    return this->GetFirstAtom()->GetRecordName();
}

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
const int& PdbResidue::GetSequenceNumber() const
{
    return this->GetFirstAtom()->GetResidueSequenceNumber();
}

const std::string& PdbResidue::GetInsertionCode() const
{
    return this->GetFirstAtom()->GetInsertionCode();
}

std::string PdbResidue::GetId() const
{
    return this->GetFirstAtom()->GetResidueId();
}

const int& PdbResidue::GetModelNumber() const
{
    return modelNumber_;
}


//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
pdb::AtomRecord* PdbResidue::FindAtom(const std::string& queryName) const
{
    for(auto &atom : atomRecords_)
    {
        if(atom->GetName() == queryName)
        {
            return atom.get();
        }
    }
    return nullptr;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidue::CreateAtom(const std::string& line, const int& currentModelNumber)
{
    atomRecords_.push_back(std::make_unique<AtomRecord>(line, currentModelNumber));
    return;
}

//void PdbResidue::AddAtom(AtomRecord* atomRecord)
//{
//    //std::cout << "Adding atom with id " << atomRecord->GetId() << " to residue.\n";
//    atomRecords_.push_back(atomRecord);
//    return;
//}

void PdbResidue::SetName(const std::string name)
{
    for(auto &atom : atomRecords_)
    {
        atom->SetResidueName(name);
    }
    return;
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidue::Print(std::ostream &out) const
{
    out << "pdb::Residue : " << this->GetId() << std::endl;
    for(auto &atom : atomRecords_)
    {
        out << "    atom : " << atom->GetId() << " X: "  << atom->GetCoordinate().GetX() << " Y: "  << atom->GetCoordinate().GetY() << " Z: " << atom->GetCoordinate().GetZ() << "\n";
    }
}

void PdbResidue::Write(std::ostream& stream) const
{
    for(auto &atom : atomRecords_)
    {
        atom->Write(stream);
    }
    if(hasTerCard_)
    {
        stream << "TER\n";
    }
}
