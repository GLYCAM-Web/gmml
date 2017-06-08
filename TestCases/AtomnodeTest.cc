#include "../includes/gmml.hpp"

using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace PdbFileSpace;

using namespace std;
int main(int argc, char* argv[]){

    //creating atomenode object
    AtomNode myatomnode;

     //creating atom pointer
     Atom *tempAtom=new Atom();

     //Initialzing Atom
     tempAtom->SetId("C1");
     tempAtom->SetElementSymbol("Chi");
     tempAtom->SetAtomType("Protein");

   //Initializng AtomNode object
    myatomnode.SetId(101);
    myatomnode.SetElementLabel("Element1");
    myatomnode.SetAtom(tempAtom);

    cout<<"Printing myatomnode object"<<endl;
    myatomnode.Print();


    cout<<"Deep Copy Operation : AtomNode myatomnode2 = myatomnode"<<endl;
    AtomNode myatomnode2 = myatomnode;

    cout<<"Printing myatomnode2 object"<<endl;
    myatomnode2.Print();

    myatomnode.GetAtom()->SetId("D1");
    cout<<"Printing myatomnode object after  myatomnode.GetAtom()->SetId(D1)"<<endl;
    myatomnode.Print();
    cout<<"Printing myatomnode2 object after  myatomnode.GetAtom()->SetId(D1)"<<endl;
    myatomnode2.Print();


    myatomnode2.GetAtom()->SetId("H1");
    cout<<"Printing myatomnode2 object after  myatomnode2.GetAtom()->SetId(H1)"<<endl;
    myatomnode2.Print();
    cout<<"Printing myatomnode object after  myatomnode2.GetAtom()->SetId(H1)"<<endl;
    myatomnode.Print();



    myatomnode.SetElementLabel("Element2");
    cout<<"Printing myatomnode object after  myatomnode.SetElementLabel(Element2)"<<endl;
    myatomnode.Print();
    cout<<"Printing myatomnode2 object after  myatomnode.SetElementLabel(Element2)"<<endl;
    myatomnode2.Print();


    myatomnode2.SetElementLabel("Element3");
    cout<<"Printing myatomnode2 object after  myatomnode.SetElementLabel(Element3)"<<endl;
    myatomnode2.Print();
    cout<<"Printing myatomnode object after  myatomnode.SetElementLabel(Element3)"<<endl;
    myatomnode.Print();


}


