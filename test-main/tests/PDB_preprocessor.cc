#include "../../includes/gmml.hpp"
#include <string>


int main(void) {

    //************************************************//
    // Details for loading in a PDB file              //
    //************************************************//

    std::vector<std::string> amino_libs, glycam_libs, other_libs, prep;
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    amino_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");

    glycam_libs.push_back("../dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib");
    glycam_libs.push_back("../dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminoct_06j_12SB.lib");
    glycam_libs.push_back("../dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib");

    other_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("../dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("../dat/CurrentParams/other/solvents.lib");

    prep.push_back("../dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j-1.prep");

    std::string parameter_file_path = "../dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j.dat";
    std::string ion_parameter_file_path = "../dat/CurrentParams/other/atomic_ions.lib";

    std::string pdb_file_path = "tests/inputs/preprocessor_input.pdb";
    PdbFileSpace::PdbFile pdbFile(pdb_file_path);

    PdbPreprocessorSpace::PdbPreprocessor preprocessor;
    preprocessor.Preprocess(&pdbFile, amino_libs, glycam_libs, other_libs, prep);
    preprocessor.ApplyPreprocessingWithTheGivenModelNumber(&pdbFile, amino_libs, glycam_libs, prep);
    pdbFile.WriteWithTheGivenModelNumber("Processed.pdb");

    return 0;
}

