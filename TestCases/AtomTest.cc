#include "../includes/gmml.hpp"

using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace PdbFileSpace;

using namespace std;
int main(int argc, char* argv[]){
    GeometryTopology::Coordinate* c_ND2 = new GeometryTopology::Coordinate(-0.847,   0.445,  -2.872);
    GeometryTopology::Coordinate* c_CG = new GeometryTopology::Coordinate(0.535,   0.914,  -3.092);


        //creating an atom
        Atom atomND2;

        atomND2.SetName("ND2");
        atomND2.AddCoordinate(c_ND2);
        atomND2.SetDescription("Atom;");
        atomND2.SetId("ND2_1_NLN_?_1_?_?_1");

        //Printing atomND2
        cout<<"Printing atomND2"<<endl;
        atomND2.Print();

     //Deep Copy without Pointer
     Atom atomCG = atomND2;

    //Deep Copy using Pointer
   //   Atom* atomCG = new Atom(*atomND2);

     //Printing atomCG
     cout<<"Printing atomCG Deep Copied"<<endl;
     atomCG.Print();

     //Changing Value of atomND2 with its Coordinate and name
     atomND2.SetName("ND3");
     atomND2.AddCoordinate(c_CG);
     cout<<"Printing atomND2 after changing Coordinate and Name"<<endl;
     atomND2.Print();

     cout<<"printing atomCG after change of atomND2"<<endl;
      atomCG.Print();

      //Changing Value of atomCG with its Coordinate and name
      atomCG.SetName("CG1");
      atomCG.AddCoordinate(c_CG);
      cout<<"Printing atomCG after changing Coordinate and Name"<<endl;
        atomCG.Print();

      cout<<"printing atomND2 after change of atomCG"<<endl;
       atomND2.Print();

   return 0;
}
