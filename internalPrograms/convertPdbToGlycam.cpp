#include "includes/gmml.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>

int main(int argc, char* argv[])
{
    // Switch for debugging statements
    // Currently positive is on, negative is off
    int local_debug = -1;

    // Usage statement and exit if incorrect number of arguments are given
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
    std::string prep = "../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep";

    MolecularModeling::Assembly assemblyA (argv[1], gmml::InputFileType::PDB);
    
    if(local_debug > 0)
    {
        // std::cout << "BuildStructureByDistance()" << std::endl;
        gmml:log(__LINE__, __FILE__, gmml::INF, "BuildStructureByDistance()");
    }
    assemblyA.BuildStructureByDistance();

    std::vector<Glycan::Monosaccharide*> monos= std::vector<Glycan::Monosaccharide*>();
    
    if(local_debug > 0)
    {
        // std::cout << "ExtractSugars()\n";
        gmml::log(__LINE__, __FILE__, gmml::INF, "ExtractSugars()");
    }
    std::vector<Glycan::Oligosaccharide*> oligos = assemblyA.ExtractSugars(amino_libs,monos,false,false);

    // Detect all side group atoms and add them to each monosaccharide
    int monoNum = 0;
    for (std::vector<Glycan::Monosaccharide*>::iterator mono_it = monos.begin(); mono_it != monos.end(); mono_it++)
    {
        if(local_debug > 0)
        {
            // std::cout << "InitiateDetectionOfCompleteSideGroupAtoms() for Mono " << monoNum << std::endl;
            gmml::log(__LINE__, __FILE__, gmml::INF, "InitiateDetectionOfCompleteSideGroupAtoms() for Mono " + std::to_string(monoNum));
            monoNum++;
        }
        (*mono_it)->InitiateDetectionOfCompleteSideGroupAtoms ();
    }

    // Get all atoms in the assembly
    AtomVector all_atoms = assemblyA.GetAllAtomsOfAssembly();

    // Update monosaccharides to make sure they are individual residues
    if(local_debug > 0)
    {
        // std::cout << "UpdateMonosaccharides2Residues()" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "UpdateMonosaccharides2Residues()");
    }
    assemblyA.UpdateMonosaccharides2Residues(monos);


    // Get the Glycam naming map
    // Initialize empty variables for ExtractResidueGlycamNamingMap()
    std::map<Glycan::Oligosaccharide*, std::vector<std::string> > oligo_id_map;
    std::map<Glycan::Oligosaccharide*, std::vector<MolecularModeling::Residue*> > oligo_residue_map;
    for (std::vector<Glycan::Oligosaccharide*>::iterator oligo_it = oligos.begin(); oligo_it != oligos.end(); oligo_it++)
    {
        std::vector<std::string> empty_vector = std::vector<std::string>();
        oligo_id_map[*oligo_it] = empty_vector;
        std::vector<MolecularModeling::Residue*> empty_residue_vector = std::vector<MolecularModeling::Residue*>();
        oligo_residue_map[*oligo_it] = empty_residue_vector;
    }

    // Get naming map and rearrange with aglycone added
    if(local_debug > 0)
    {
        // std::cout << "ExtractResidueGlycamNamingMap()" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "ExtractResidueGlycamNamingMap()");
    }
    gmml::GlycamResidueNamingMap res_map = assemblyA.ExtractResidueGlycamNamingMap(oligos, oligo_id_map, oligo_residue_map);

    if(local_debug > 0)
    {
        // std::cout << "PutAglyconeInNewResidueAndRearrangeGlycanResidues)" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "PutAglyconeInNewResidueAndRearrangeGlycanResidues()");
    }
    assemblyA.PutAglyconeInNewResidueAndRearrangeGlycanResidues(oligos, oligo_residue_map);

    if(local_debug > 0)
    {
        // std::cout << "TestUpdateResidueName2GlycamName()" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "TestUpdateResidueName2GlycamName()");
    }
    assemblyA.TestUpdateResidueName2GlycamName(res_map, prep);

    //Match and rename atoms
    //std::map<MolecularModeling::Atom*, MolecularModeling::Atom*> actual_template_match;
    std::map<Glycan::Oligosaccharide*, pdb2glycam_matching_tracker*> match_tracker;
    if(local_debug > 0)
    {
        // std::cout << "MatchPdbAtoms2Glycam()" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "MatchPdbAtoms2Glycam()");
    }
    assemblyA.MatchPdbAtoms2Glycam(oligo_residue_map, prep, match_tracker);
    if(local_debug > 0)
    {
        // std::cout << "Done with MatchPdbAtoms2Glycam()" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Done with MatchPdbAtoms2Glycam()");
    }

    // Assign name from naming map to oligosaccharides
    for (unsigned int i = 0; i < oligos.size(); i++)
    {
        if(local_debug > 0)
        {
            // std::cout << "Assigning Name to Oligosaccharide: " << i << std::endl;
            gmml::log(__LINE__, __FILE__, gmml::INF, "Assigning Name to Oligosaccharide: " + std::to_string(i));
        }
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
            for (std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>::iterator mapit = first_match.begin(); mapit != first_match.end(); mapit++)
            {
                mapit->first->SetName(mapit->second->GetName());
            }
        }
    }
    if(local_debug > 0)
    {
        // std::cout << "BuildPdbFileStructureFromAssembly()" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "BuildPdbFileStructureFromAssembly()");
    }
    PdbFileSpace::PdbFile *outputPdbFile = assemblyA.BuildPdbFileStructureFromAssembly();
    if (argc == 3)
    {
        if(local_debug > 0)
        {
            // std::cout << "outputPdbFile->Write()" << std::endl;
            gmml::log(__LINE__, __FILE__, gmml::INF, "outputPdbFile->Write()");
        }
        std::stringstream outNameStream;
        outNameStream << argv[2] << ".pdb";
        outputPdbFile->Write(outNameStream.str());
    }
    else
    {
        if(local_debug > 0)
        {
            // std::cout << "outputPdbFile->Write()" << std::endl;
            gmml::log(__LINE__, __FILE__, gmml::INF, "outputPdbFile->Write()");
        }
        outputPdbFile->Write("pdb2glycam_output.pdb");
    }

    if(local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Done with pdb2glycam!");
    }
}
