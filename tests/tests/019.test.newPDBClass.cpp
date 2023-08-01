#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include <string>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << " tests/inputs/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    pdb::PdbFile pdbFile(argv[1]);
    pdb::PreprocessorOptions options; // Default values are good.
    std::cout << "Preprocessing\n";
    pdb::PreprocessorInformation ppInfo = pdbFile.PreProcess(options);
    for (auto& assembly : pdbFile.getAssemblies()) // Just testing, doing it this way to get around const in Ensemble.
                                                   // ToDo: Why is there a const blockage in Ensemble?
    {
        std::cout << "Bonding atoms by distance for assembly" << std::endl;
        ;
        cds::bondAtomsByDistance(assembly->getAtoms());
        // OFF molecule
        try
        {
            std::ofstream outFileStream;
            outFileStream.open("outputOffFile.off");
            std::string outFileName = "Assembly";
            cds::WriteAssemblyToOffFile(assembly, outFileStream, outFileName.c_str());
            outFileStream.close();
        }
        catch (std::runtime_error& error)
        {
            std::stringstream ss;
            ss << "Runtime error thrown when writing to off file:\n" << error.what() << "\n";
            std::cout << ss.str();
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        }
        catch (...)
        {
            std::cout << "Unknown error when writing to off file.\n";
            gmml::log(__LINE__, __FILE__, gmml::ERR, "Unknown error when writing to off file.\n");
        }
    }
    std::cout << "Finished bonding atoms by distance" << std::endl;
    pdbFile.Write("./outputPdbFile.pdb");

    // Just showing what's in the ppInfo and how to access it
    std::cout << "Unrecognized atoms:\n";
    for (auto& unrecognized : ppInfo.unrecognizedAtoms_)
    {
        std::cout << unrecognized.name_ << " | " << unrecognized.residue_.getName() << " | "
                  << unrecognized.residue_.getChainId() << " | " << unrecognized.residue_.getNumberAndInsertionCode()
                  << "\n";
    }
    std::cout << "Missing heavy atoms:\n";
    for (auto& missing : ppInfo.missingHeavyAtoms_)
    {
        std::cout << missing.name_ << " | " << missing.residue_.getName() << " | " << missing.residue_.getChainId()
                  << " | " << missing.residue_.getNumberAndInsertionCode() << "\n";
    }
    std::cout << "Unrecognized residues:\n";
    for (auto& unrecognized : ppInfo.unrecognizedResidues_)
    {
        std::cout << unrecognized.getName() << " | " << unrecognized.getChainId() << " | "
                  << unrecognized.getNumberAndInsertionCode() << "\n";
    }
    std::cout << "Gaps in amino acid chain:\n";
    for (auto& gap : ppInfo.missingResidues_)
    {
        std::cout << gap.chainId_ << " | " << gap.residueBeforeGap_ << " | " << gap.residueAfterGap_ << " | "
                  << gap.terminationBeforeGap_ << " | " << gap.terminationAfterGap_ << "\n";
    }
    std::cout << "Histidine Protonation:\n";
    for (auto& his : ppInfo.hisResidues_)
    {
        std::cout << his.getName() << " | " << his.getChainId() << " | " << his.getNumberAndInsertionCode() << "\n";
    }
    std::cout << "Disulphide bonds:\n";
    for (auto& cysBond : ppInfo.cysBondResidues_)
    {
        std::cout << cysBond.residue1_.getChainId() << " | " << cysBond.residue1_.getName() << " | "
                  << cysBond.residue1_.getNumberAndInsertionCode() << " | " << cysBond.distance_ << " | "
                  << cysBond.residue2_.getChainId() << " | " << cysBond.residue2_.getName() << " | "
                  << cysBond.residue2_.getNumberAndInsertionCode() << "\n";
    }
    std::cout << "Chain terminations:\n";
    for (auto& chainT : ppInfo.chainTerminals_)
    {
        std::cout << chainT.chainId_ << " | " << chainT.startIndex_ << " | " << chainT.nTermination_ << " | "
                  << chainT.endIndex_ << " | " << chainT.cTermination_ << "\n";
    }
    std::cout << "NonNatural Protein Residues:\n";
    for (auto& nonNaturalResidue : ppInfo.nonNaturalProteinResidues_)
    {
        std::cout << nonNaturalResidue.residue_.getChainId() << " | " << nonNaturalResidue.residue_.getName() << " | "
                  << nonNaturalResidue.residue_.getNumberAndInsertionCode() << "\n";
    }
    return 0;
}
