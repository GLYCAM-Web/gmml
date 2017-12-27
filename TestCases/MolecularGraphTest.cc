#include "../includes/gmml.hpp"

using namespace MolecularModeling;
using namespace GeometryTopology;


using namespace std;
int main(int argc, char* argv[]){

   //creating five atoms object of Atom class
    Atom a1;
    Atom *a1Pointer= &a1;

    Atom a2;
    Atom *a2Pointer= &a2;

    Atom a3;
    Atom *a3Pointer= &a3;

    Atom a4;
    Atom *a4Pointer= &a4;

    Atom a5;
    Atom *a5Pointer= &a5;



// Creating node of type atom
    Node node1;
    node1.SetNode(a1Pointer);
    node1.SetNodeId(1);
    Node *node1ptr=&node1;

    Node node2;
    node2.SetNode(a2Pointer);
    node2.SetNodeId(2);
    Node *node2ptr=&node2;


    Node node3;
    node3.SetNode(a3Pointer);
    node3.SetNodeId(3);
    Node *node3ptr=&node3;


    Node node4;
    node4.SetNode(a4Pointer);
    node4.SetNodeId(4);
    Node *node4ptr=&node4;

    Node node5;
    node5.SetNode(a5Pointer);
    node5.SetNodeId(5);
    Node *node5ptr=&node5;


//Creating Graph
    MolecularGraph mygraph;


    mygraph.AddNewNode(node1ptr);
    mygraph.AddNewNode(node2ptr);
    mygraph.AddNewNode(node3ptr);
    mygraph.AddNewNode(node4ptr);
    mygraph.AddNewNode(node5ptr);



mygraph.AddEdge(node1ptr,node5ptr);
mygraph.AddEdge(node5ptr,node4ptr);
mygraph.AddEdge(node4ptr,node2ptr);
mygraph.AddEdge(node2ptr,node3ptr);


mygraph.Print();


return 0;
}

