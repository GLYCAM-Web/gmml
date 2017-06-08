#include "../includes/gmml.hpp"

using namespace MolecularModeling;
using namespace std;


int main(int argc, char* argv[]){

   cout<<"*******Deep Copy Without Pointer*******"<<endl;
   cout<<endl;
   Element e1;
   e1.SetSymbol("E1");
   e1.SetName("Element E1");

   //Deep Copy
   Element e2=e1;

    cout<<"Printing element E1"<<endl;
   e1.Print();
    cout<<endl;

   cout<<"Printing element E2 after Deep Copy"<<endl;
   e2.Print();
    cout<<endl;
    //Changing Name of Element E1  to E3
    e1.SetName("Element E3");

    cout<<"Printing element E1 after name change to E3"<<endl;
    e1.Print();
     cout<<endl;

   cout<<"Printing element E2 after name change of E1 to E3. E2 still hold the previous name as E1"<<endl;
   e2.Print();
 cout<<endl;

   //Changing Name of Element E2 to E4
    e2.SetName("Element E4");

    cout<<"Printing element E2 after name change to E4"<<endl;
    e2.Print();
    cout<<endl;

    cout<<"Printing element E1 after name change of E2 to E4.E1 still hold the previous name as E3"<<endl;
    e1.Print();


    //Deep Copy using Pointer
     //Element *e3=new Element(*e1);
}
