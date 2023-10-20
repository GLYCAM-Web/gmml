#include <string>
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"

#include <fstream> // std::ifstream

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << " tests/inputs/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    // requirement: a chain ID for every single ATOM entry, and all ligand atoms should be put in a single residue.
    pdb::PdbFile pdbFile(argv[1]); // PdbFile is an "Ensemble" (made up of "Assemblies"), but if you want to just set
                                   // every molecule to have any chain ID you can do:
    for (auto& residue : pdbFile.getResidues())
    {
        static_cast<pdb::PdbResidue*>(residue)->setChainId("Y");
        // std::cout << "Set chain of " << residue->getStringId() << "\n";
    }
    // ResidueTypes are guessed upon input. Using that guess to find the ligand, can improve this if you need:
    std::vector<cds::Residue*> ligandResidues =
        cdsSelections::selectResiduesByType(pdbFile.getResidues(), cds::ResidueType::Undefined);
    if (ligandResidues.empty())
    {
        std::cout << "No ligand residues found in input file\n";
        return 0;
    }
    cds::Residue* firstLigandResidue = ligandResidues.front();
    for (auto& ligandResidue : ligandResidues) // Each MODEL in PdbFile is converted into an "Assembly"
    {                                          // Every ligand residue gets the same residue number as the first one.
        // std::cout << "Renumbering and renaming " << ligandResidue->getStringId() << "\n";
        ligandResidue->setNumber(firstLigandResidue->getNumber());
        ligandResidue->setName(firstLigandResidue->getName());
    }
    pdbFile.Write("./026.outputPdbFile.pdb");
    // Separate thing showing this functionality to Yao.
    // Don't want the selection updated every frame so do it outside the loop:#

    cds::Ensemble selectionEnsemble;
    for (auto& assembly : pdbFile.getAssemblies())
    {
        std::vector<cds::Residue*> myResidues = assembly->getResidues();
        cds::Residue* queryResidue            = codeUtils::findElementWithNumber(
            myResidues, 5); // somehow you specify this in inputs. e.g. A_405 chain A, residue 405.
        double distanceFromQueryResidue             = 10.3; // inputs
        std::vector<cds::Residue*> selectedResidues = cdsSelections::selectResiduesWithinDistanceN(
            myResidues, queryResidue, distanceFromQueryResidue); // Need to write function.
        cds::Assembly newAssembly(selectedResidues);             // Need to add constructor
        selectionEnsemble.addAssembly(newAssembly);              // done
    }
    const std::string outName = "026.outputSelection.pdb";
    std::ofstream outFileStream;
    outFileStream.open(outName.c_str());
    cds::writeEnsembleToPdb(outFileStream, selectionEnsemble.getAssemblies());
    outFileStream.close();
    pdbFile.Write("./026.outputPdbFileWithoutSelection.pdb");
    return 0;
}
