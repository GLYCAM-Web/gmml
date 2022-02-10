#include "includes/gmml.hpp"
#include "includes/InputSet/PdbFile/pdbFile.hpp"
#include "includes/Resolver/NewPdbPreprocessor/pdbPreprocessorInputs.hpp"

#include <string>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << "tests/inputs/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    pdb::PdbFile pdbFile(argv[1]);
    pdb::PreprocessorOptions options; // Default values are good.
    pdb::PreprocessorInformation ppInfo = pdbFile.PreProcess(options);
    pdbFile.Write("./outputPdbFile.pdb");
    // Just showing what's in the ppInfo and how to access it
    std::cout << "Unrecognized atoms:\n";
    for(auto &unrecognized : ppInfo.unrecognizedAtoms_)
    {
        std::cout << unrecognized.name_ << " | " << unrecognized.residue_.name_ << " | " << unrecognized.residue_.chainId_ << " | " << unrecognized.residue_.number_ << "\n";
    }
    std::cout << "Missing heavy atoms:\n";
    for(auto &missing : ppInfo.missingHeavyAtoms_)
    {
        std::cout << missing.name_ << " | " << missing.residue_.name_ << " | " << missing.residue_.chainId_ << " | " << missing.residue_.number_ << "\n";
    }
    std::cout << "Unrecognized residues:\n";
    for(auto &unrecognized : ppInfo.unrecognizedResidues_)
    {
        std::cout << unrecognized.name_ << " | " << unrecognized.chainId_ << " | " << unrecognized.number_ << "\n";
    }
    std::cout << "Gaps in amino acid chain:\n";
    for(auto &gap : ppInfo.missingResidues_)
    {
        std::cout << gap.chainId_ << " | " << gap.residueBeforeGap_ << " | " << gap.residueAfterGap_ << " | " << gap.terminationBeforeGap_ << " | " << gap.terminationAfterGap_ << "\n";
    }
    std::cout << "Histidine Protonation:\n";
    for(auto &his : ppInfo.hisResidues_)
    {
        std::cout << his.name_ << " | " << his.chainId_ << " | " << his.number_ << "\n";
    }
    std::cout << "Disulphide bonds:\n";
    for(auto &cysBond : ppInfo.cysBondResidues_)
    {
        std::cout << cysBond.residue1_.chainId_ << " | " <<  cysBond.residue1_.name_ << " | " <<  cysBond.residue1_.number_ << " | " <<  cysBond.distance_ << " | " << cysBond.residue2_.chainId_ << " | " <<  cysBond.residue2_.name_ << " | " <<  cysBond.residue2_.number_<< "\n";
    }
    std::cout << "Chain terminations:\n";
    for(auto &chainT : ppInfo.chainTerminals_)
    {
        std::cout << chainT.chainId_ << " | " << chainT.startIndex_ << " | " << chainT.nTermination_ <<  " | " << chainT.endIndex_ << " | " << chainT.cTermination_ << "\n";
    }
    return 0;
}
