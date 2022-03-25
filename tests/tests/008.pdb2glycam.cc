#include "includes/gmml.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>

int main(int argc, char* argv[])
{
    if ( (argc != 2) && (argc != 3) )
    {
        std::cout << "Usage: pdb2glycam inputFile.pdb [outputFileName]\n";
        std::cout << "Example: pdb2glycam 1RVX.pdb 1RVX_GLYCAMIFICATION.pdb\n";
        std::cout << "outFileName is optional, the default name is pdb2glycam_output.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    typedef std::vector<Glycan::Oligosaccharide*> OligosaccharideVector;
    std::vector<std::string> amino_libs;
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");
    std::string prep = "../dat/prep/GLYCAM_06j-1.prep";
    MolecularModeling::Assembly assemblyA (argv[1], gmml::InputFileType::PDB);
    //std::cout << "BuildStructureByDistance()" << std::endl;
    assemblyA.BuildStructureByDistance();
    std::vector<Glycan::Monosaccharide*> monos= std::vector<Glycan::Monosaccharide*>();
    //std::cout << "ExtractSugars\n";
    std::vector<Glycan::Oligosaccharide*> oligos = assemblyA.ExtractSugars(amino_libs,monos,false,false);
    for (std::vector<Glycan::Monosaccharide*>::iterator mono_it = monos.begin(); mono_it != monos.end(); mono_it++){
        //std::cout << "InitiateDetectionOfCompleteSideGroupAtoms()" << std::endl;
        (*mono_it)->InitiateDetectionOfCompleteSideGroupAtoms ();
    }
    AtomVector all_atoms = assemblyA.GetAllAtomsOfAssembly();
    //std::cout << "UpdateMonosaccharides2Residues()" << std::endl;
    assemblyA.UpdateMonosaccharides2Residues(monos);
    std::map<Glycan::Oligosaccharide*, std::vector<std::string> > oligo_id_map;
    std::map<Glycan::Oligosaccharide*, std::vector<MolecularModeling::Residue*> > oligo_residue_map;
    for (std::vector<Glycan::Oligosaccharide*>::iterator oligo_it = oligos.begin(); oligo_it != oligos.end(); oligo_it++){
        std::vector<std::string> empty_vector = std::vector<std::string>();
        oligo_id_map[*oligo_it] = empty_vector;
        std::vector<MolecularModeling::Residue*> empty_residue_vector = std::vector<MolecularModeling::Residue*>();
        oligo_residue_map[*oligo_it] = empty_residue_vector;
    }
    //std::cout << "ExtractResidueGlycamNamingMap()" << std::endl;
    gmml::GlycamResidueNamingMap res_map = assemblyA.ExtractResidueGlycamNamingMap(oligos, oligo_id_map, oligo_residue_map);
    assemblyA.PutAglyconeInNewResidueAndRearrangeGlycanResidues(oligos);
    //std::cout << "TestUpdateResidueName2GlycamName()" << std::endl;
    assemblyA.TestUpdateResidueName2GlycamName(res_map, prep);
    //Match and rename atoms
    std::map<MolecularModeling::Atom*, MolecularModeling::Atom*> actual_template_match;
    assemblyA.MatchPdbAtoms2Glycam(oligo_residue_map, prep, actual_template_match);
    for (std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>::iterator mapit = actual_template_match.begin(); mapit != actual_template_match.end(); mapit++){
        mapit->first->SetName(mapit->second->GetName());
    }
    //std::cout << "BuildPdbFileStructureFromAssembly()" << std::endl;
    PdbFileSpace::PdbFile *outputPdbFile = assemblyA.BuildPdbFileStructureFromAssembly();
    if (argc == 3)
    {
        std::stringstream outNameStream;
        outNameStream << argv[2] << ".pdb";
        outputPdbFile->Write(outNameStream.str());
    }
    else
    {
        //std::cout << "outputPdbFile->Write()" << std::endl;
        outputPdbFile->Write("pdb2glycam_output.pdb");
    }
}
