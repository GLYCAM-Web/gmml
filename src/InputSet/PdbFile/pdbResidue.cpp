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
pdb::AtomRecord* PdbResidue::GetFirstAtom() const
{
    if (atomRecords_.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "atomRecords in PdbResidue is empty. Impossible!?");
        std::exit(1);
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

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
pdb::AtomRecord* PdbResidue::FindAtom(const std::string selector) const
{
    pdb::AtomRecord* nullAtom = nullptr;
    for(auto &atom : atomRecords_)
    {
        std::size_t found = atom->GetId().find(selector);
        if (found != std::string::npos)
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
    std::cout << "Adding atom with id " << atomRecord->GetId() << " to residue.\n";
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
