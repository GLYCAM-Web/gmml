#include "../includes/gmml.hpp"
using namespace MolecularModeling;

using namespace std;
int main(int argc, char* argv[]){

    MolecularDynamicAtom matom1;
    matom1.SetAtomType("hydrogen");
    matom1.Print();

    cout<<"Deep Copying:  MolecularDynamicAtom matom2=matom1;"<<endl;
    MolecularDynamicAtom matom2=matom1;

    matom2.Print();

      cout<<"Changing matom1 Atom Type form hydrogen to iron"<<endl;
     matom1.SetAtomType("iron");
     matom1.Print();


    cout<<"matom2 Atom Type still holds previous value as hydrogen"<<endl;
    matom2.Print();


    cout<<"Changing matom2 Atom Type form hydrogen to chlorine"<<endl;
    matom2.SetAtomType("chlorine");
     matom2.Print();

     cout<<"matom1 Atom Type still holds previous value as iron"<<endl;
     matom1.Print();
}
