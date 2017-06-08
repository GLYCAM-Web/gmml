#include "includes/gmml.hpp"
#include <string>
#include <iostream>

using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace PdbFileSpace;

using namespace std;
int main(int argc, char* argv[]){
    GeometryTopology::Coordinate* c_ND2 = new GeometryTopology::Coordinate(-0.847,   0.445,  -2.872);
    GeometryTopology::Coordinate* c_CG = new GeometryTopology::Coordinate(0.535,   0.914,  -3.092);


   Atom test2;
   test2.SetName("Ayu");
test2.AddCoordinate(c_CG);
test2.Print();
   Atom atomz = test2;
    test2.SetName("Aj");
    cout<<endl;
   atomz.Print();
      
   atomz.SetName("jais");
atomz.AddCoordinate(c_ND2);
 atomz.Print();
test2.Print();
   return 0;
}
