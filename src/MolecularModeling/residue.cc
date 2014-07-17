#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/atom.hpp"

using namespace std;
using namespace MolecularModeling;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Residue::Residue() {}

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
    atoms_ = atoms;
}
void Residue::AddAtom(Atom *atom)
{
    atoms_.push_back(atom);
}
void Residue::SetHeadAtoms(AtomVector head_atoms)
{
    head_atoms_ = head_atoms;
}
void Residue::AddHeadAtom(Atom *head_atom)
{
    head_atoms_.push_back(head_atom);
}
void Residue::SetTailAtoms(AtomVector tail_atoms)
{
    tail_atoms_ = tail_atoms;
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

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Residue::Print(ostream &out)
{
}
