#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/CodeUtils/logging.hpp"
using pdb::PdbResidue;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidue::PdbResidue(AtomRecord* atomRecord)
{
    this->AddAtom(atomRecord);
}
PdbResidue::PdbResidue(std::vector<AtomRecord*> atomRecords)
{
    atomRecords_ = atomRecords;
}
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
    return atomRecords_.back();
}

pdb::AtomRecord* PdbResidue::GetFirstAtom() const
{
    if (atomRecords_.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "atomRecords in PdbResidue is empty. Impossible!?");
    }
    return atomRecords_.front();
}
const std::string& PdbResidue::GetChainId() const
{
    return this->GetFirstAtom()->GetChainId();
}
const std::string& PdbResidue::GetName() const
{
    return this->GetFirstAtom()->GetResidueName();
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
    return this->GetFirstAtom()->GetModelNumber();
}


//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
pdb::AtomRecord* PdbResidue::FindAtom(const std::string& queryName) const
{
    for(auto &atom : atomRecords_)
    {
//        std::size_t found = atom->GetId().find(name);
//        if (found != std::string::npos)
//        {
//            return atom;
//        }
        if(atom->GetName() == queryName)
        {
            return atom;
        }
    }
    return nullptr;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidue::AddAtom(AtomRecord* atomRecord)
{
    //std::cout << "Adding atom with id " << atomRecord->GetId() << " to residue.\n";
    atomRecords_.push_back(atomRecord);
    return;
}
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
