#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/residueproperties.hpp"
#include "../../includes/MolecularModeling/residuenode.hpp"
#include "../../../includes/MolecularModeling/overlaps.hpp"
#include "../../includes/common.hpp"
#include <algorithm>    // std::any_of

//using namespace std; // It is forbidden.
using namespace MolecularModeling;
using namespace gmml;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Residue::Residue() {}

Residue::Residue(Assembly *assembly, std::string name)
{
    assembly_ = assembly;
    name_ = name;
    atoms_ = AtomVector();
    head_atoms_ = AtomVector();
    tail_atoms_ = AtomVector();
    chemical_type_ = "";
    description_ = "";
    id_ = "";
}

Residue::Residue(Residue *residue)
{
    assembly_ = new Assembly(residue->GetAssembly());
    name_ = residue->GetName();
    atoms_ = AtomVector();
    AtomVector atoms = residue->GetAtoms();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
        atoms_.push_back(new Atom(*it));
    head_atoms_ = AtomVector();
    AtomVector head_atoms = residue->GetHeadAtoms();
    for(AtomVector::iterator it = head_atoms.begin(); it != head_atoms.end(); it++)
        head_atoms_.push_back(new Atom(*it));
    tail_atoms_ = AtomVector();
    AtomVector tail_atoms = residue->GetTailAtoms();
    for(AtomVector::iterator it = tail_atoms.begin(); it != tail_atoms.end(); it++)
        tail_atoms_.push_back(new Atom(*it));
    chemical_type_ = residue->GetChemicalType();
    description_ = residue->GetDescription();
    id_ = residue->GetId();
}

Residue::Residue(Residue& residue){

    Assembly *tempAssembly=residue.GetAssembly();
    this->assembly_=tempAssembly;

    this->name_=residue.GetName();
    this->atoms_ = AtomVector();

    AtomVector atoms = residue.GetAtoms();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
        atoms_.push_back(new Atom(*it));

    this->head_atoms_ = AtomVector();
    AtomVector head_atoms = residue.GetHeadAtoms();
    for(AtomVector::iterator it = head_atoms.begin(); it != head_atoms.end(); it++)
        head_atoms_.push_back(new Atom(*it));

    this->tail_atoms_ = AtomVector();
    AtomVector tail_atoms = residue.GetTailAtoms();
    for(AtomVector::iterator it = tail_atoms.begin(); it != tail_atoms.end(); it++)
        tail_atoms_.push_back(new Atom(*it));

    this->chemical_type_ = residue.GetChemicalType();
    this->description_ = residue.GetDescription();
    this->id_ = residue.GetId();

}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Assembly* Residue::GetAssembly()
{
    return assembly_;
}
std::string Residue::GetName()
{
    return name_;
}
std::string Residue::GetNumber()
{
    StringVector id = gmml::Split(id_, "_");
    return id.at(2); // This is silly, why not add residue number to class? OG: I know right?
}
Residue::AtomVector Residue::GetAtoms()
{
    return atoms_;
}
Residue::AtomVector Residue::GetHeadAtoms()
{
    return head_atoms_;
}
Residue::AtomVector Residue::GetTailAtoms()
{
    return tail_atoms_;
}
std::string Residue::GetChemicalType()
{
    return chemical_type_;
}
std::string Residue::GetDescription()
{
    return description_;
}
std::string Residue::GetId()
{
    return id_;
}

 //Added by ayush on 11/20/17 for residuenode in assembly
ResidueNode* Residue::GetNode()
{
        return node_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Residue::SetAssembly(Assembly *assembly)
{
    assembly_ = assembly;
}
void Residue::SetName(std::string name)
{
    name_ = name;
}
void Residue::SetAtoms(AtomVector atoms)
{
    atoms_.clear();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_.push_back(*it);
    }
}
void Residue::AddAtom(Atom *atom)
{
    atoms_.push_back(atom);
}
void Residue::RemoveAtom(Atom *atom)
{
    AtomVector newAtoms = AtomVector();
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* a = *it;
        if(a->GetId().compare(atom->GetId()) != 0)
        {
            if(a->GetNode() != NULL)
                a->GetNode()->RemoveNodeNeighbor(atom);
            newAtoms.push_back(a);
        }
    }
    this->SetAtoms(newAtoms);
}

void Residue::SetHeadAtoms(AtomVector head_atoms)
{
    head_atoms_.clear();
    for(AtomVector::iterator it = head_atoms.begin(); it != head_atoms.end(); it++)
    {
        head_atoms_.push_back(*it);
    }
}
void Residue::AddHeadAtom(Atom *head_atom)
{
    head_atoms_.push_back(head_atom);
}
void Residue::SetTailAtoms(AtomVector tail_atoms)
{
    tail_atoms_.clear();
    for(AtomVector::iterator it = tail_atoms.begin(); it != tail_atoms.end(); it++)
    {
        tail_atoms_.push_back(*it);
    }
}
void Residue::AddTailAtom(Atom *tail_atom)
{
    tail_atoms_.push_back(tail_atom);
}
void Residue::SetChemicalType(std::string chemical_type)
{
    chemical_type_ = chemical_type;
}
void Residue::SetDescription(std::string description)
{
    description_ = description;
}
void Residue::SetId(std::string id)
{
    id_ = id;
}
void Residue::ReplaceAtomCoordinates(AtomVector *newAtoms)
{
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); ++it)
    {
        Atom *atom = (*it);
        for(AtomVector::iterator itt = newAtoms->begin(); itt != newAtoms->end(); ++itt)
        {
            Atom *atom1 = (*itt);
            //std::cout << "Comparing with " << atom1->GetName() << std::endl;
            if (atom->GetName() == atom1->GetName() )
            {
                //std::cout << "Replacing " << atom1->GetName() << " with " << atom->GetName() << std::endl;
                //std::cout << "Before X=" << atom->GetCoordinates().at(0)->GetX() << std::endl;
                atom->GetCoordinates().at(0)->SetX( atom1->GetCoordinates().at(0)->GetX() );
                atom->GetCoordinates().at(0)->SetY( atom1->GetCoordinates().at(0)->GetY() );
                atom->GetCoordinates().at(0)->SetZ( atom1->GetCoordinates().at(0)->GetZ() );
                //std::cout << "After X=" << atom->GetCoordinates().at(0)->GetX() << std::endl;
            }
        }
    }
}

 //Added by ayush on 11/20/17 for residuenode in assembly
void Residue::SetNode(ResidueNode* node)
{
    node_ = node;
}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
bool Residue::CheckSymbolBasedElementLabeling()
{
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        if(atom->GetElementSymbol().compare("") == 0)
            return false;
    }
    return true;
}

bool Residue::CheckParameterBasedElementLabeling()
{
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        if(atom->GetAtomType().compare("UNK") == 0 || atom->GetAtomType().compare("") == 0)
            return false;
    }
    return true;
}

bool Residue::GraphElementLabeling()
{
    if(this->CheckSymbolBasedElementLabeling())
        return this->GraphSymbolBasedElementLabeling();
    else if(this->CheckParameterBasedElementLabeling())
        return this->GraphParameterBasedElementLabeling();
    return GraphPredictionBasedElementLabeling();
}

bool Residue::GraphSymbolBasedElementLabeling()
{
    std::cout << "Labeling residue nodes based on elements' symbol ... " << std::endl;
    bool flag = true;
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
            atom_node->SetElementLabel(atom->GetElementSymbol());
        else
        {
            atom_node->SetElementLabel("UNK");
            flag = false;
        }
    }
    return flag;
}

bool Residue::GraphParameterBasedElementLabeling()
{
    std::cout << "Labeling residue nodes based on atom type and parameter file" << std::endl;
    bool flag = true;
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            std::string element_symbol = gmml::AtomTypesLookup(atom->GetAtomType()).element_symbol_;
            if(element_symbol != "")
                atom_node->SetElementLabel(element_symbol);
            else
            {
                atom_node->SetElementLabel("UNK");
                flag = false;
            }
        }
        else
        {
            atom_node->SetElementLabel("UNK");
            flag = false;
        }
    }
    return flag;
}

bool Residue::GraphPredictionBasedElementLabeling()
{
    std::cout << "Labeling residue nodes based on first letter prediction ... " << std::endl;
    bool flag = true;
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
            atom_node->SetElementLabel(isdigit(atom->GetName().at(0)) ? atom->GetName().at(1) + "" : atom->GetName().at(0) + "");
        else
        {
            atom_node->SetElementLabel("UNK");
            flag = false;
        }
    }
    return flag;
}

Residue::AtomVector Residue::GetAtomsWithLowestIntraDegree()
{
    int degree = INFINITY;
    AtomVector lowest_degree_atoms = AtomVector();
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        int current = atom->GetNode()->GetIntraEdgeDegree();
        if(current < degree)
            current = degree;
    }
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        if(atom->GetNode()->GetIntraEdgeDegree() == degree)
            lowest_degree_atoms.push_back(atom);
    }

    return lowest_degree_atoms;
}

double Residue::CalculateAtomicOverlaps(Assembly *assemblyB)
{
    AtomVector assemblyBAtoms = assemblyB->GetAllAtomsOfAssembly();
    AtomVector residueAtoms = this->GetAtoms();
    return gmml::CalculateAtomicOverlaps(residueAtoms, assemblyBAtoms);
}
double Residue::CalculateAtomicOverlaps(AtomVector assemblyBAtoms)
{
    AtomVector residueAtoms = this->GetAtoms();
    return gmml::CalculateAtomicOverlaps(residueAtoms, assemblyBAtoms);
}

bool Residue::CheckIfProtein()
{
    if( std::find( PROTEINS, ( PROTEINS + PROTEINSSIZE ), this->GetName() ) != ( PROTEINS + PROTEINSSIZE ) )
    {
        return true;
    }
    return false;
}

bool Residue::CheckIfWater() {
	if( this->GetName().compare( "HOH" ) == 0 ) {
		return true;
	}
	return false;
}

/*GeometryTopology::Coordinate Residue::GetRingCenter() // Disabled by OG; GetIsRing returns true for all atoms even when IsRing wasn't set. 
{
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    int numberOfRingAtoms = 0;
    AtomVector atoms = this->GetAtoms();

    for(Assembly::AtomVector::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
    {
        if ( (*atom)->GetIsRing() )
        {
            numberOfRingAtoms++;
            //std::cout << "Atom is ring: " << (*atom)->GetName() << std::endl;
            sumX += (*atom)->GetCoordinates().at(0)->GetX();
            sumY += (*atom)->GetCoordinates().at(0)->GetY();
            sumZ += (*atom)->GetCoordinates().at(0)->GetZ();
        }
    }
    GeometryTopology::Coordinate center;
    center.SetX( sumX / numberOfRingAtoms  );
    center.SetY( sumY / numberOfRingAtoms  );
    center.SetZ( sumZ / numberOfRingAtoms  );
    return center;
}*/

GeometryTopology::Coordinate Residue::GetGeometricCenter()
{
    AtomVector atoms = this->GetAtoms();
    if(atoms.size() == 0)
    {
        std::cout << "Problem in Residue::GetGeometricCenter(), the residue " << this->GetId() << " contains no atoms." << std::endl;
    }

    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    for(Assembly::AtomVector::iterator atom = atoms.begin(); atom != atoms.end(); ++atom)
    {   
       // std::cout << "atoms size is " << atoms.size() << " for " << (*atom)->GetId() << std::endl;
        sumX += (*atom)->GetCoordinates().at(0)->GetX();
        sumY += (*atom)->GetCoordinates().at(0)->GetY();
        sumZ += (*atom)->GetCoordinates().at(0)->GetZ();
    }
    GeometryTopology::Coordinate center;
    center.SetX( sumX / atoms.size()  );
    center.SetY( sumY / atoms.size()  );
    center.SetZ( sumZ / atoms.size()  );
    return center;
}

Atom* Residue::GetAtom(std::string query_name)
{
    Atom* return_atom;
    AtomVector atoms = this->GetAtoms();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
        if ((*it)->GetName().compare(query_name)==0)
        {
            return_atom = (*it);
        }
    }
    return return_atom; // may be unset
}

Atom* Residue::GetAtom(unsigned long long query_index)
{
    Atom* return_atom;
    AtomVector atoms = this->GetAtoms();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
        if ((*it)->GetIndex() == query_index)
        {
            return_atom = (*it);
        }
    }
    return return_atom; // may be unset
}

Atom* Residue::GetAtom(int query_id)
{
    Atom* return_atom;
    AtomVector atoms = this->GetAtoms();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
        if ((*it)->GetId() == query_id)
        {
            return_atom = (*it);
        }
    }
    return return_atom; // may be unset
}



//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Residue::Print(std::ostream &out)
{
    out << "------------------------ " << name_ << " --------------------------" << std::endl;
    out << "Head atoms: ";
    for(AtomVector::iterator it = head_atoms_.begin(); it != head_atoms_.end(); it++)
    {
        Atom* atom = *it;
//        out << atom->GetResidue()->GetName() << ":" << atom->GetName() << "; ";
        out << atom->GetId() << "; ";
    }
    out << std::endl;
    out << "Tail atoms: ";
    for(AtomVector::iterator it = tail_atoms_.begin(); it != tail_atoms_.end(); it++)
    {
        Atom* atom = *it;
//        out << atom->GetResidue()->GetName() << ":" << atom->GetName() << "; ";
        out << atom->GetId() << "; ";
    }
    out << std::endl;
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        atom->Print(out);
    }
}

void Residue::PrettyPrintHet(std::ostream &out)
{
    out << "------------------------ " << "Residue " << " --------------------------" << std::endl;
    out << " ID: " << id_ << std::endl;
    out << " Name: " << name_ << std::endl;
//    out << " Chemical type: " << chemical_type_ << std::endl;
//    out << " Description: " << description_ << std::endl;
    out << " ATOMS: ";

    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        out << atom->GetId() << ", ";
    }

    out << std::endl;
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        out << "------------------------ " << "Atom" << " --------------------------" << std::endl;
        out << " ID: " << atom->GetId() << std::endl;
        out << " Name: " << atom->GetName() << std::endl;
        out << " Atom type: " << atom->GetAtomType() << std::endl;
        out << " Charge: " << atom->GetCharge() << std::endl;
//        out << " Chemical Type: " << atom->GetChemicalType() << std::endl;
//        out << " Description: " << atom->GetDescription() << std::endl;
        out << " Mass: " << atom->GetMass() << std::endl;
        GeometryTopology::Coordinate* coords = atom->GetCoordinates().at(0);
        out << " Coordinates" <<  " X: " << coords->GetX() << ", Y: " << coords->GetY() << ", Z: " << coords->GetZ() << std::endl;
        out << " Neighbors: ";
        AtomNode* node = atom->GetNode();
        AtomVector neighbors = node->GetNodeNeighbors();
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            Atom* neighbor = *it1;
            out << neighbor->GetId() << ", ";
        }
        out << std::endl;
        atom->PrintHet(out);
    }
}

void Residue::PrintHetResidues(std::ostream &out)
{
    out << id_ << ";" << name_ << ";";
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        out << atom->GetId() << ",";
    }
    out << std::endl;
}

void Residue::PrintHetAtoms(std::ostream &out)
{
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        out << atom->GetId() << ";" << atom->GetName() << ";" << atom->GetAtomType() << ";" << atom->GetCharge() << ";" << atom->GetMass() << ";";
        GeometryTopology::Coordinate* coords = atom->GetCoordinates().at(0);
        out << coords->GetX() << "," << coords->GetY() << "," << coords->GetZ() << ";";

        AtomNode* node = atom->GetNode();
        AtomVector neighbors = node->GetNodeNeighbors();
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            Atom* neighbor = *it1;
            out << neighbor->GetId() << ",";
        }
        out << std::endl;
    }
}

void Residue::WriteHetResidues(std::ofstream& out)
{
    out << id_ << ";" << name_ << ";";
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        if (it == (atoms_.end() - 1) )
            out << atom->GetId();
        else
            out << atom->GetId() << ",";
    }
    out << std::endl;
}

void Residue::WriteHetAtoms(std::ofstream& out)
{
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        out << atom->GetId() << ";" << atom->GetName() << ";" << atom->GetAtomType() << ";" << atom->GetCharge() << ";" << atom->GetMass() << ";";
        GeometryTopology::Coordinate* coords = atom->GetCoordinates().at(0);
        out << coords->GetX() << "," << coords->GetY() << "," << coords->GetZ() << ";";

        AtomNode* node = atom->GetNode();
        AtomVector neighbors = node->GetNodeNeighbors();
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            Atom* neighbor = *it1;
            if(it1 == (neighbors.end() - 1))
                out << neighbor->GetId();
            else
                out << neighbor->GetId() << ",";
        }
        out << std::endl;
    }
}
