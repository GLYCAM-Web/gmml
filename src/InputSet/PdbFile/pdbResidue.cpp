#include <includes/InputSet/PdbFile/pdbAtom.hpp>
#include <sstream>
#include <string>
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
#include "includes/InputSet/PdbFile/pdbResidueId.hpp" // residueId
#include "includes/GeometryTopology/geometrytopology.hpp" // get_cartesian_point_from_internal_coords
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp" //RemoveWhiteSpace

using pdb::PdbResidue;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
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
            this->addAtom(std::make_unique<PdbAtom>(line));
        }
    }
    this->SetType( this->determineType( this->getName() ) );
    return;
}

PdbResidue::PdbResidue(const std::string residueName, const PdbResidue *referenceResidue)
: cds::Residue(residueName, referenceResidue)
{ // should instead call copy constructor and then rename with residueName?
//    this->setName(residueName); // handled by cdsResidue cTor
//    this->setNumber(referenceResidue->getNumber()); // handled by cdsResidue cTor
    this->setInsertionCode(referenceResidue->getInsertionCode());
    this->setChainId(referenceResidue->getChainId());
    this->SetType(referenceResidue->GetType());
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

const std::string PdbResidue::getNumberAndInsertionCode() const
{
    return std::to_string(this->getNumber()) + this->getInsertionCode();
}

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
        const cds::Atom* atom = this->FindAtom("H");
        if (atom != nullptr)
        {
            gmml::log(__LINE__,__FILE__,gmml::INF, "Deleting atom with id: " + static_cast<const PdbAtom*>(atom)->GetId());
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
        const cds::Atom* atom = this->FindAtom("OXT");
        if (atom == nullptr)
        {
            // I don't like this, but at least it's somewhat contained:
            const cds::Atom* atomCA = this->FindAtom("CA");
            const cds::Atom* atomC = this->FindAtom("C");
            const cds::Atom* atomO = this->FindAtom("O");
            GeometryTopology::Coordinate oxtCoord = GeometryTopology::get_cartesian_point_from_internal_coords(atomCA->getCoordinate(), atomC->getCoordinate(), atomO->getCoordinate(), 120.0, 180.0, 1.25);
            this->addAtom(std::make_unique<PdbAtom>("OXT", oxtCoord));
            gmml::log(__LINE__,__FILE__,gmml::INF, "Created new atom named OXT after " + static_cast<const PdbAtom*>(atomO)->GetId());
        }
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    }
    return;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidue::Print(std::ostream &out) const
{
    out << "pdb::Residue : " << this->printId() << std::endl;
    for(auto &atom : this->getAtoms())
    {
        out << "    atom : " << static_cast<const PdbAtom*>(atom)->GetId() << " X: "  << atom->getCoordinate()->GetX() << " Y: "  << atom->getCoordinate()->GetY() << " Z: " << atom->getCoordinate()->GetZ() << "\n";
    }
}

void PdbResidue::Write(std::ostream& stream) const
{
    this->WritePdb(stream, this->HasTerCard());
}
