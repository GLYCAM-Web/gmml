#include "../includes/gmml.hpp"
#include "../includes/MolecularModeling/atom.hpp"
#include "../includes/MolecularModeling/Molecule.hpp"
#include "../includes/MolecularModeling/residuenode.hpp"

using namespace MolecularModeling;
using namespace GeometryTopology;
using namespace PdbFileSpace;
using namespace TopologyFileSpace;

using namespace std;
int main(int argc, char* argv[]){

    cout<<"This is Topology File operation!!"<<endl;
    std::vector<std::string> amino_libs, glycam_libs, other_libs, prep;
    amino_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    amino_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    amino_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");
    glycam_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib");
    glycam_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminoct_06j_12SB.lib");
    glycam_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib");
    other_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("/home/ayush/gems/gmml/dat/CurrentParams/other/solvents.lib");
    prep.push_back("/home/ayush/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j-1.prep");
    std::string parameter_file_path = "/home/ayush/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j.dat";
    std::string ion_parameter_file_path = "/home/ayush/gems/gmml/dat/CurrentParams/other/atomic_ions.lib";

    std::string sequence = "DNeup5Aca2-6DGalpNAca1-OME";
   // std::string outputFileName;
    cout<<sequence<<endl;
    cout<< prep.at(0)<<endl;
    cout<<parameter_file_path<<endl;
    Assembly assembly;
    //CondensedSequence *condensedSequence = new CondensedSequence(sequence);
    //CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_info = condensedSequence->GetCondensedSequenceRotamersAndGlycosidicAnglesInfo(condensedSequence->GetCondensedSequenceResidueTree());
    cout<<"Before -->BuildAssemblyFromCondensedSequence"<<endl;
    assembly.BuildAssemblyFromCondensedSequence(sequence, prep.at(0),parameter_file_path, true);

    cout<<"After -->BuildAssemblyFromCondensedSequence"<<endl;
//cout<<"Before -->AddSolvent"<<endl;
  //  assembly->AddSolvent(5, 1,assembly, "/programs/amber16/dat/leap/lib/tip3pbox.off");
    //cout<<"After -->AddSolvent"<<endl;

   typedef std::vector<ResidueNode*> ResidueNodeVector;
    ResidueNodeVector allAsseblyResidueNode = assembly.GenerateResidueNodesInAssembly();

         for(ResidueNodeVector::iterator it3 = allAsseblyResidueNode.begin(); it3 != allAsseblyResidueNode.end(); it3++)
         {
                ResidueNode* currentResidueNode = (*it3);
                currentResidueNode->Print();
         }
         assembly.GenerateMoleculesInAssembly();

//     for(AtomVector::iterator it3 = allAsseblyAtom.begin(); it3 != allAsseblyAtom.end(); it3++)
//     {
//         Atom* currentAtom = *it3;

//         Residue r1;

//         cout<<"Head Atom "<<currentAtom->GetIndex()<<" Residue it belong: "<<currentAtom->GetResidue()->GetName()<<endl;

//         AtomNode* currentAtomNode = currentAtom->GetNode();
//         AtomVector currentAtomNodeVector= currentAtomNode->GetNodeNeighbors();
//         if(!currentAtomNodeVector.empty())
//         {
//                cout<<"Neighbour Atoms of head"<<endl;
//             for(AtomVector::iterator it4 = currentAtomNodeVector.begin(); it4 != currentAtomNodeVector.end(); it4++)
//             {
//                    Atom* currentAtom = *it4;
//                 cout<<" "<<currentAtom->GetIndex()<<" Residue it belong: "<<currentAtom->GetResidue()->GetName()<<endl;
//             }
//         }else{

//             cout<<"Head Atom does not have any neighbours"<<endl;
//         }


//         cout<<endl;
//         cout<<"--------------------------------------------------------------------------------------------------------------"<<endl;


//     }

//     typedef std::vector<Molecule*> MoleculeVector;

//    MoleculeVector tempMolecules = assembly.GetMolecules();
//     for(MoleculeVector::iterator it = tempMolecules.begin(); it != tempMolecules.end(); it++)
//     {
//         Molecule* molecule = (*it);
//         cout<<"Index:"<<molecule->GetMoleculeIndex()<<" "<<molecule->GetMoleculeType()<<endl;
//     }


//    TopologyFileSpace::TopologyFile *outputTopologyFile1 = assembly->BuildTopologyFileStructureFromAssembly(parameter_file_path, ion_parameter_file_path);
//    outputTopologyFile1->Write("/home/ayush/Desktop/test-BuildFromSequence.parm7");

//    CoordinateFileSpace::CoordinateFile *outputCoordinateFile = assembly.BuildCoordinateFileStructureFromAssembly();
//    outputCoordinateFile->Write("/home/ayush/Desktop/test-BuildFromSequence.rst7");

}



