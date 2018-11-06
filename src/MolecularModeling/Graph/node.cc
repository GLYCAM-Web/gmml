#include "../../includes/MolecularModeling/Graph/node.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/MolecularModeling/molecule.hpp"


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

   GraphDS::Node::Node(): is_visited_(false){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
MolecularModeling::Atom* GraphDS::Node::GetAtomNode()
{
    return node_atom_;
}

MolecularModeling::Residue* GraphDS::Node::GetResidueNode()
{
    return node_residue_;
}

MolecularModeling::Molecule* GraphDS::Node::GetMoleculeNode()
{
    return node_molecule_;
}

std::string GraphDS::Node::GetNodeId()
{
    return node_id_;
}

bool GraphDS::Node::GetIsVisited()
{
    return is_visited_;
}

GraphDS::Node::EdgeVector GraphDS::Node::GetEdgeList()
{
    return edgeList_;
}

GraphDS::Node::TagsVector GraphDS::Node::GetTagList()
{
    return tags_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

void GraphDS::Node::SetAtomNode(MolecularModeling::Atom* node_atom)
{
    this->node_atom_= node_atom;
}

void GraphDS::Node::SetResidueNode(MolecularModeling::Residue* node_residue)
{
    this->node_residue_= node_residue;
}

void GraphDS::Node::SetMoleculeNode(MolecularModeling::Molecule *node_molecule)
{
    this->node_molecule_= node_molecule;
}

void GraphDS::Node::SetNodeId(long node_id)
{
    this->node_id_=node_id;
}

void GraphDS::Node::SetIsVisited(bool is_visited)
{
    this->is_visited_=is_visited;
}

void GraphDS::Node::SetEdgeList(EdgeVector edgeList)
{
    this->edgeList_.clear();
    for(EdgeVector::iterator it = edgeList.begin(); it != edgeList.end(); it++) {
        this->edgeList_.push_back(*it);
    }
}

void GraphDS::Node::SetTagList(TagsVector tags)
{
    this->tags_.clear();
    for(TagsVector::iterator it = tags.begin(); it != tags.end(); it++) {
        this->tags_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

void GraphDS::Node::AddEdge(Edge *edge)
{
    edgeList_.push_back(edge);
}
void GraphDS::Node::AddTag(std::string tag)
{
    tags_.push_back(tag);
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

void GraphDS::Node::Print(std::ostream &out)
{}
