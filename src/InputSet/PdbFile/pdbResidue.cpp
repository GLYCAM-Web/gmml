#include <sstream>
#include <string>
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/InputSet/PdbFile/pdbResidueId.hpp" // residueId
#include "includes/InputSet/PdbFile/pdbFunctions.hpp" // extractResidueId
#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp" // get_cartesian_point_from_internal_coords
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp" //RemoveWhiteSpace

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
//PdbResidue::PdbResidue(const std::string& line, const int& modelNumber)
//{
//    //atoms_.push_back(std::make_unique<AtomRecord>(line, modelNumber));
//    this->CreateAtomFromLine(line, modelNumber);
//    modelNumber_ = modelNumber;
//    hasTerCard_ = false;
//}
PdbResidue::PdbResidue(std::stringstream &singleResidueSecion, std::string firstLine)
{
    ResidueId resId(firstLine);
    this->setName(resId.getName());
    this->setNumber(std::stoi(resId.getNumber()));
    this->setInsertionCode(resId.getInsertionCode());
    this->setChainId(resId.getChainId());
    std::string line;
    while(getline(singleResidueSecion, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        if ( (recordName == "ATOM") || (recordName == "HETATM") )
        {
            this->addAtom(std::make_unique<AtomRecord>(line));
        }
    }
    return;
}

PdbResidue::PdbResidue(const std::string residueName, const PdbResidue *referenceResidue)
: cds::cdsResidue<AtomRecord>(residueName, referenceResidue)
{ // should instead call copy constructor and then rename with residueName?
//    this->setName(residueName); // handled by cdsResidue cTor
//    this->setNumber(referenceResidue->getNumber()); // handled by cdsResidue cTor
    this->setInsertionCode(referenceResidue->getInsertionCode());
    this->setChainId(referenceResidue->getChainId());
    return;
}

//PdbResidue::PdbResidue(std::vector<std::unique_ptr<AtomRecord>>& atomRecords)
//{
//    atomRecordss_.swap(atomRecords); // atomRecords will become empty, atoms_ will contain contents of atomRecords.
//    hasTerCard_ = false;
//}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
pdb::ResidueId PdbResidue::getId() const
{
    ResidueId temp(this->getName(), std::to_string(this->getNumber()), this->getInsertionCode(), this->getChainId());
    return temp;
}
const std::string& PdbResidue::getNumberAndInsertionCode() const
{
    return std::to_string(this->getNumber()) + this->getInsertionCode();
}
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
        return "N" + this->getName();
    }
    else if (this->containsLabel("CTerminal"))
    {
        return "C" + this->getName();
    }
    return this->getName();
}
//const int& PdbResidue::GetSequenceNumber() const
//{
//    return this->GetFirstAtom()->GetResidueSequenceNumber();
//}

//const std::string& PdbResidue::GetId() const
//{
//    return this->getName() + "_" + this->getInsertionCode() + "_" + this->getChainId() + "_" + this->getModelNumber();
//}


//const int& PdbResidue::GetModelNumber() const
//{
//    return modelNumber_;
//}
//////////////////////////////////////////////////////////
//                    MUTATOR                           //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void PdbResidue::modifyNTerminal(const std::string& type)
{
    gmml::log(__LINE__,__FILE__,gmml::INF, "Modifying N Terminal of : " + this->printId());
    if (type == "NH3+")
    {
        AtomRecord* atom = this->FindAtom("H");
        if (atom != nullptr)
        {
            gmml::log(__LINE__,__FILE__,gmml::INF, "Deleting atom with id: " + atom->GetId());
            this->deleteAtom(atom);
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}

void PdbResidue::modifyCTerminal(const std::string& type)
{
    gmml::log(__LINE__,__FILE__,gmml::INF, "Modifying C Terminal of : " + this->printId());
    if (type == "CO2-")
    {
        AtomRecord* atom = this->FindAtom("OXT");
        if (atom == nullptr)
        {
            // I don't like this, but at least it's somewhat contained:
            AtomRecord* atomCA = this->FindAtom("CA");
            AtomRecord* atomC = this->FindAtom("C");
            AtomRecord* atomO = this->FindAtom("O");
            GeometryTopology::Coordinate oxtCoord = GeometryTopology::get_cartesian_point_from_internal_coords(atomCA->GetCoordinate(), atomC->GetCoordinate(), atomO->GetCoordinate(), 120.0, 180.0, 1.25);
            this->createAtom("OXT", oxtCoord);
            gmml::log(__LINE__,__FILE__,gmml::INF, "Created new atom named OXT after " + atomO->GetId());
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}

//void PdbResidue::CreateAtomFromLine(const std::string& line, const int& currentModelNumber)
//{
//    AtomRecord tempRecord(line, currentModelNumber);
//    this->createAtom(tempRecord);
//    //atoms_.push_back(std::make_unique<AtomRecord>(line, currentModelNumber));
//    return;
//}


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
    out << "pdb::Residue : " << this->printId() << std::endl;
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
