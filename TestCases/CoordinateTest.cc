#include "../includes/gmml.hpp"

using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace PdbFileSpace;

using namespace std;
int main(int argc, char* argv[]){

// Deep Copy using Pointer

cout<<"Performing Deep Copy using Pointer"<<endl;
Coordinate *c1=new Coordinate();
c1->SetX(1.25);
c1->SetY(1.50);
c1->SetZ(1.75);
cout<<"Printing Coordinate c1"<<endl;
c1->Print();
cout<<endl;


Coordinate *c2=new Coordinate(c1);
cout<<"Printing Coordinate c2 Deep Copied from c1"<<endl;
c2->Print();
cout<<endl;
cout<<"Changing x values of Coordinate c1"<<endl;
c1->SetX(2.222);
cout<<"Printing Coordinate c1"<<endl;
c1->Print();
cout<<endl;
cout<<"Printing Coordinate c2"<<endl;
c2->Print();
cout<<endl;

cout<<"Changing x values of Coordinate c2"<<endl;
c2->SetX(4.444);
cout<<"Printing Coordinate c2"<<endl;
c2->Print();
cout<<endl;
cout<<"Printing Coordinate c1"<<endl;
c1->Print();
cout<<endl;
cout<<"----------------------------------------------------------------------"<<endl;
cout<<" "<<endl;

cout<<"Performing Deep Copy without Pointer"<<endl;
Coordinate c3;
c3.SetX(5.989);
c3.SetY(3.65);
c3.SetZ(4.647);
cout<<"Printing Coordinate c3"<<endl;
c3.Print();
cout<<endl;

Coordinate c4=c3;
cout<<"Printing Coordinate c2 Deep Copied from c1"<<endl;
c4.Print();
cout<<endl;

cout<<"Changing x values of Coordinate c3"<<endl;
c3.SetX(8.888);
cout<<"Printing Coordinate c3"<<endl;
c3.Print();
cout<<endl;

cout<<"Printing Coordinate c4"<<endl;
c4.Print();
cout<<endl;

cout<<"Changing x values of Coordinate c4"<<endl;
c4.SetX(1.1111);
cout<<"Printing Coordinate c3"<<endl;
c3.Print();
cout<<endl;

cout<<"Printing Coordinate c4"<<endl;
c4.Print();
cout<<endl;


}


