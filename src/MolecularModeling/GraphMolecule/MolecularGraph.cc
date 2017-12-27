#include "../../../includes/MolecularModeling/GraphMolecule/Node.hpp"
#include "../../../includes/MolecularModeling/GraphMolecule/MolecularGraph.hpp"


using namespace std;
using namespace MolecularModeling;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
MolecularGraph::MolecularGraph(){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void MolecularGraph::AddNewNode(Node *newNode)
{
    graphNodeList_.push_back(newNode);
}


void MolecularGraph::RemoveNode(Node* delNode)
{
    for(NodeVector::iterator it = graphNodeList_.begin(); it != graphNodeList_.end(); it++)
    {
        Node* node = (*it);
        if(node->GetNodeId()==(delNode->GetNodeId()))
        graphNodeList_.erase(it);
    }
}


void MolecularGraph::AddEdge(Node* firstNode, Node* secondNode)
{
    firstNode->AddadjNode(secondNode);
    secondNode->AddadjNode(firstNode);
}

Node* MolecularGraph::FindNodeById(long node_id)
{

}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void MolecularGraph::Print(ostream &out)
{
     out << "------------------------ Displaying The Graph : --------------------------" << endl;

     for(NodeVector::iterator it = graphNodeList_.begin(); it != graphNodeList_.end(); it++)
     {
            Node *node= *it;
            out<<"Head:"<<node->GetNodeId()<<"-> ";
            node->Print(out);
            out<<endl;
     }

}

