#include "../includes/gmml.hpp"

using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace PdbFileSpace;

using namespace std;
int main(int argc, char* argv[]){

Residue r;
//r.Print();
r.SetChemicalType("chi");
r.SetName("Ayush");
Residue r2=r;

r2.Print();
cout<<"r2name"<<r2.GetName()<<endl;
r2.SetName("Jaiswal");
cout<<"r2name"<<r2.GetName()<<endl;
}
