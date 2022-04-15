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
    typedef MolecularModeling::Assembly::Pdb2glycamMatchingTracker pdb2glycam_matching_tracker;
    typedef MolecularModeling::Assembly::Pdb2glycamMatchingFailInfo pdb2glycam_matching_fail_info;

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
    assemblyA.PutAglyconeInNewResidueAndRearrangeGlycanResidues(oligos, oligo_residue_map);
    //std::cout << "TestUpdateResidueName2GlycamName()" << std::endl;
    assemblyA.TestUpdateResidueName2GlycamName(res_map, prep);
    //Match and rename atoms
    //std::map<MolecularModeling::Atom*, MolecularModeling::Atom*> actual_template_match;
    std::map<Glycan::Oligosaccharide*, pdb2glycam_matching_tracker*> match_tracker;
    assemblyA.MatchPdbAtoms2Glycam(oligo_residue_map, prep, match_tracker);

    for (unsigned int i = 0; i < oligos.size(); i++){
        pdb2glycam_matching_tracker* this_oligo_match_tracker = match_tracker[oligos[i]];
        std::vector<std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>>& all_isomorphisms = this_oligo_match_tracker->all_isomorphisms;

        if (all_isomorphisms.empty()){
            std::cout << "Oligosaccharide " << i+1 << " matching failed." << std::endl;
            std::cout << "Here are the atoms that might be the issue:" << std::endl;

            int largest_iteration_length = this_oligo_match_tracker->largest_iteration_length;
            std::cout << "Largest iteration length: " << largest_iteration_length << std::endl;

            std::vector<pdb2glycam_matching_fail_info*>& failures = this_oligo_match_tracker->failures;
            for (unsigned int j = 0; j < failures.size(); j++){
                pdb2glycam_matching_fail_info* this_failure = failures[j];
                if (this_failure->iteration_length == largest_iteration_length){
                    std::cout << "Failed atom: " << this_failure->failed_atom->GetResidue()->GetName() << "-" << this_failure->failed_atom->GetName() << std::endl;
                    std::cout << "Failure message: " << this_failure->failure_notice << std::endl << std::endl;
                }
            }

        }
        else{
            std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>& first_match = all_isomorphisms[0];
            for (std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>::iterator mapit = first_match.begin(); mapit != first_match.end(); mapit++){
                mapit->first->SetName(mapit->second->GetName());
            }
        }
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
