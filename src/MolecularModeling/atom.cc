#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/quantommechanicatom.hpp"
#include "../../includes/MolecularModeling/moleculardynamicatom.hpp"
#include "../../includes/MolecularModeling/dockingatom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/residue.hpp"

using namespace std;
using namespace MolecularModeling;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Atom::Atom() : name_(""), chemical_type_(""), element_symbol_(""), description_(""), id_("")
{
    coordinates_ = CoordinateVector();
    residue_ = NULL;
    node_ = NULL;
}

Atom::Atom(Residue *residue, string name, CoordinateVector coordinates) :
    chemical_type_(""), element_symbol_(""), description_("")
{
    residue_ = residue;
    name_ = name;
    coordinates_ = CoordinateVector();
    for(CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
        coordinates_.push_back(*it);
    node_ = NULL;
}

Atom::Atom(Atom *atom)
{
    residue_ = new Residue(atom->GetResidue());
    name_ = atom->GetName();
    coordinates_ = CoordinateVector();
    CoordinateVector coordinates = atom->GetCoordinates();
    for(CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
        coordinates_.push_back(new GeometryTopology::Coordinate(*it));

    AtomNode node = atom->GetNode();
    node_ = new AtomNode(node);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
Residue* Atom::GetResidue()
{
    return residue_;
}
string Atom::GetName()
{
    return name_;
}
Atom::CoordinateVector Atom::GetCoordinates()
{
    return coordinates_;
}
string Atom::GetChemicalType()
{
    return chemical_type_;
}
string Atom::GetDescription()
{
    return description_;
}
string Atom::GetElementSymbol()
{
    return element_symbol_;
}
AtomNode* Atom::GetNode()
{
    return node_;
}
string Atom::GetId()
{
    return id_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Atom::SetResidue(Residue *residue)
{
    residue_ = residue;
}
void Atom::SetName(string name)
{
    name_ = name;
}
void Atom::SetCoordinates(CoordinateVector coordinates)
{
    coordinates_.clear();
    for(CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
    {
        coordinates_.push_back(*it);
    }
}
void Atom::AddCoordinate(GeometryTopology::Coordinate *coordinate)
{
    coordinates_.push_back(coordinate);
}
void Atom::SetChemicalType(string chemical_type)
{
    chemical_type_ = chemical_type;
}
void Atom::SetDescription(string description)
{
    description_ = description;
}
void Atom::SetElementSymbol(string element_symbol)
{
    element_symbol_ = element_symbol;
}
void Atom::SetNode(AtomNode *node)
{
    node_ = node;
}
void Atom::SetId(string id)
{
    id_ = id;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Atom::Print(ostream &out)
{
    out << "Atom name: " << name_ << endl;
    out << "Coordinates: " << endl;
    if(coordinates_.size() != 0)
    {
        for(CoordinateVector::iterator it = coordinates_.begin(); it != coordinates_.end(); it++)
        {
            GeometryTopology::Coordinate* coordinate = *it;
            out << "\t";
            coordinate->Print(out);
            out << endl;
        }
    }
    out << "**************** Structure *****************" << endl;
    if(node_ != NULL)
        node_->Print(out);
    out << endl;

}
