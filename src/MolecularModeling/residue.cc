#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"

using namespace std;
using namespace MolecularModeling;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Residue::Residue() {}

Residue::Residue(Assembly *assembly, string name)
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

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Assembly* Residue::GetAssembly()
{
    return assembly_;
}
string Residue::GetName()
{
    return name_;
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
string Residue::GetChemicalType()
{
    return chemical_type_;
}
string Residue::GetDescription()
{
    return description_;
}
string Residue::GetId()
{
    return id_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Residue::SetAssembly(Assembly *assembly)
{
    assembly_ = assembly;
}
void Residue::SetName(string name)
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
void Residue::SetChemicalType(string chemical_type)
{
    chemical_type_ = chemical_type;
}
void Residue::SetDescription(string description)
{
    description_ = description;
}
void Residue::SetId(string id)
{
    id_ = id;
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
    cout << "Labeling residue nodes based on elements' symbol ... " << endl;
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
    cout << "Labeling residue nodes based on atom type and parameter file" << endl;
    bool flag = true;
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        AtomNode* atom_node = atom->GetNode();
        if(atom_node != NULL)
        {
            string element_symbol = gmml::AtomTypesLookup(atom->GetAtomType()).element_symbol_;
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
    cout << "Labeling residue nodes based on first letter prediction ... " << endl;
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

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Residue::Print(ostream &out)
{
    out << "------------------------ " << name_ << " --------------------------" << endl;
    out << "Head atoms: ";
    for(AtomVector::iterator it = head_atoms_.begin(); it != head_atoms_.end(); it++)
    {
        Atom* atom = *it;
//        out << atom->GetResidue()->GetName() << ":" << atom->GetName() << "; ";
        out << atom->GetId() << "; ";
    }
    out << endl;
    out << "Tail atoms: ";
    for(AtomVector::iterator it = tail_atoms_.begin(); it != tail_atoms_.end(); it++)
    {
        Atom* atom = *it;
//        out << atom->GetResidue()->GetName() << ":" << atom->GetName() << "; ";
        out << atom->GetId() << "; ";
    }
    out << endl;
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        atom->Print(out);
    }
}

void Residue::PrettyPrintHet(ostream &out)
{
    out << "------------------------ " << "Residue " << " --------------------------" << endl;
    out << " ID: " << id_ << endl;
    out << " Name: " << name_ << endl;
//    out << " Chemical type: " << chemical_type_ << endl;
//    out << " Description: " << description_ << endl;
    out << " ATOMS: ";

    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        out << atom->GetId() << ", ";
    }

    out << endl;
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        out << "------------------------ " << "Atom" << " --------------------------" << endl;
        out << " ID: " << atom->GetId() << endl;
        out << " Name: " << atom->GetName() << endl;
        out << " Atom type: " << atom->GetAtomType() << endl;
        out << " Charge: " << atom->GetCharge() << endl;
//        out << " Chemical Type: " << atom->GetChemicalType() << endl;
//        out << " Description: " << atom->GetDescription() << endl;
        out << " Mass: " << atom->GetMass() << endl;
        GeometryTopology::Coordinate* coords = atom->GetCoordinates().at(0);
        out << " Coordinates" <<  " X: " << coords->GetX() << ", Y: " << coords->GetY() << ", Z: " << coords->GetZ() << endl;
        out << " Neighbors: ";
        AtomNode* node = atom->GetNode();
        AtomVector neighbors = node->GetNodeNeighbors();
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            Atom* neighbor = *it1;
            out << neighbor->GetId() << ", ";
        }
        out << endl;
        atom->PrintHet(out);
    }
}

void Residue::PrintHetResidues(ostream &out)
{
    out << id_ << ";" << name_ << ";";
    for(AtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        Atom* atom = *it;
        out << atom->GetId() << ",";
    }
    out << endl;
}

void Residue::PrintHetAtoms(ostream &out)
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
        out << endl;
    }
}

void Residue::WriteHetResidues(ofstream& out)
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
    out << endl;
}

void Residue::WriteHetAtoms(ofstream& out)
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
        out << endl;
    }
}
