#include "../includes/gmml.hpp"

using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace PdbFileSpace;

using namespace std;
int main(int argc, char* argv[]){

    //Min Coordinate
    Coordinate *min=new Coordinate();
    min->SetX(1.25);
    min->SetY(1.50);
    min->SetZ(1.75);

    //Max Coordinate
    Coordinate *max=new Coordinate();
    max->SetX(9.25);
    max->SetY(9.50);
    max->SetZ(9.75);

    Grid grid1;
    Assembly *assembly1 =new Assembly();
    assembly1->SetChemicalType("Iron");
    grid1.SetAssembly(assembly1);

    Grid::CellVector grid_cells = Grid::CellVector();

  grid1.SetCells(grid_cells);

    grid1.SetMinCorner(min);
    grid1.SetMaxCorner(max);
   // grid1.Print();

    grid1.UpdateGrid(2,3,4);

   Grid grid2=grid1;

//    cout<<grid2.GetMinCorner()->GetX()<<endl;
//   grid2.GetMinCorner()->SetX(3.33);

//   cout<<grid2.GetMinCorner()->GetX()<<endl;

//   cout<<"grid1.minx"<<grid1.GetMinCorner()->GetX()<<endl;

//   grid2.GetMinCorner()->SetX(1.33);
//   cout<<"grid1.minx"<<grid2.GetMinCorner()->GetX()<<endl;
//   cout<<"grid1.minx"<<grid1.GetMinCorner()->GetX()<<endl;

//    Assembly *assembly1 ;
//    assembly1->SetChemicalType("Iron");
//    grid1.SetAssembly(assembly1);

cout<<"grid1 "<<grid1.GetAssembly()->GetChemicalType()<<endl;
cout<<"grid2 "<<grid2.GetAssembly()->GetChemicalType()<<endl;

grid2.GetAssembly()->SetChemicalType("Protein");
cout<<"grid2 "<<grid2.GetAssembly()->GetChemicalType()<<endl;
cout<<"grid1 "<<grid1.GetAssembly()->GetChemicalType()<<endl;


grid1.GetAssembly()->SetChemicalType("Glucose");
cout<<"grid1 "<<grid1.GetAssembly()->GetChemicalType()<<endl;
cout<<"grid2 "<<grid2.GetAssembly()->GetChemicalType()<<endl;

//cout<<"Assembly1-grid1 :"<<grid1.GetAssembly()->GetChemicalType()<<endl;

//cout<<"Assembly1-grid2 :"<<grid2.GetAssembly()->GetChemicalType()<<endl;
//

}


