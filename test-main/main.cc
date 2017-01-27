#include "includes/gmml.hpp"
#include <string>

using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace PdbFileSpace;


int main(int argc, char* argv[]){
    GeometryTopology::Coordinate* c_ND2 = new GeometryTopology::Coordinate(-0.847,   0.445,  -2.872);
    GeometryTopology::Coordinate* c_CG = new GeometryTopology::Coordinate(0.535,   0.914,  -3.092);

    Atom *atomND2 = new Atom();
    atomND2->SetName("ND2");
    atomND2->AddCoordinate(c_ND2);
    atomND2->SetDescription("Atom;");
    atomND2->SetId("ND2_1_NLN_?_1_?_?_1");

    Atom* atomCG = new Atom();
    atomCG->SetName("CG");
    atomCG->AddCoordinate(c_CG);
    atomCG->SetDescription("Atom;");
    atomCG->SetId("CG_2_NLN_?_1_?_?_1");

    Residue* residue_NLN = new Residue();
    residue_NLN->SetName("NLN");
    residue_NLN->AddAtom(atomCG);
    residue_NLN->AddAtom(atomND2);
    residue_NLN->SetId("NLN_?_1_?_?_1");

    Assembly* assembly_NLN = new Assembly();
    assembly_NLN->AddResidue(residue_NLN);

    atomND2->SetResidue(residue_NLN);
    atomCG->SetResidue(residue_NLN);

    residue_NLN->SetAssembly(assembly_NLN);

    assembly_NLN->BuildStructureByDistance();
    assembly_NLN->Print();
    PdbFileSpace::PdbFile *outputPdbFileNLN = assembly_NLN->BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFileNLN->Write("test-NLN.pdb");
    return 0;
}
