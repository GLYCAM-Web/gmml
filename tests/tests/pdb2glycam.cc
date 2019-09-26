#include "../../includes/gmml.hpp"
#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "../../includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "../../includes/ParameterSet/OffFileSpace/offfileprocessingexception.hpp"
#include "../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../includes/InputSet/Utilities/response.hpp"
#include "../../includes/Glycan/note.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>

int main(int argc, char* argv[])
{
    typedef std::vector<Glycan::Oligosaccharide*> OligosaccharideVector;
    std::vector<std::string> amino_libs;
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");
    std::string prep = "../dat/prep/GLYCAM_06j-1.prep";
    MolecularModeling::Assembly assemblyA (argv[1], gmml::InputFileType::PDB);

    assemblyA.BuildStructureByDistance();
    std::vector<Glycan::Monosaccharide*> monos= std::vector<Glycan::Monosaccharide*>();
    std::vector<Glycan::Oligosaccharide*> oligos = assemblyA.ExtractSugars(amino_libs,monos,false,false);

    for (std::vector<Glycan::Monosaccharide*>::iterator mono_it = monos.begin(); mono_it != monos.end(); mono_it++){
        (*mono_it)->InitiateDetectionOfCompleteSideGroupAtoms ();
    }

    AtomVector all_atoms = assemblyA.GetAllAtomsOfAssembly();
    assemblyA.UpdateMonosaccharides2Residues(monos);

    std::map<Glycan::Oligosaccharide*, std::vector<std::string> > oligo_id_map;
    std::map<Glycan::Oligosaccharide*, std::vector<MolecularModeling::Residue*> > oligo_residue_map;

    for (std::vector<Glycan::Oligosaccharide*>::iterator oligo_it = oligos.begin(); oligo_it != oligos.end(); oligo_it++){
        std::vector<std::string> empty_vector = std::vector<std::string>();
        oligo_id_map[*oligo_it] = empty_vector;
        std::vector<MolecularModeling::Residue*> empty_residue_vector = std::vector<MolecularModeling::Residue*>();
        oligo_residue_map[*oligo_it] = empty_residue_vector;
    }

    gmml::GlycamResidueNamingMap res_map = assemblyA.ExtractResidueGlycamNamingMap(oligos, oligo_id_map, oligo_residue_map);
    assemblyA.TestUpdateResidueName2GlycamName(res_map, prep);
    //assemblyA.UpdateResidueName2GlycamName(res_map, prep);

    assemblyA.RenameAtoms(oligo_residue_map, prep);
    PdbFileSpace::PdbFile *outputPdbFile = assemblyA.BuildPdbFileStructureFromAssembly();
    outputPdbFile->Write("pdb2glycam_output.pdb");
}
