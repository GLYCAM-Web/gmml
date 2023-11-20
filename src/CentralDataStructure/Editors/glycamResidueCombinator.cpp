#include "includes/CentralDataStructure/Editors/glycamResidueCombinator.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06LinkageCodes.hpp"

std::vector<std::string> residueCombinator::selectAllAtomsThatCanBeSubstituted(std::vector<cds::Atom*> queryAtoms)
{
    std::vector<std::string> foundNames;
    for (auto& atom : queryAtoms)
    {
        if (atom->getType() == "Oh") // If a hydroxyl.
        {
            foundNames.push_back(atom->getName().substr(1, 1));
            std::cout << atom->getName() << ", ";
        }
    }
    std::cout << std::endl;
    return foundNames;
}

// ChatGPT made this:
std::vector<std::vector<std::string>> residueCombinator::getCombinations(const std::vector<std::string>& elements)
{
    std::vector<std::vector<std::string>> combinations;
    std::vector<std::string> currentCombination;
    // Recursive function to generate all combinations
    std::function<void(size_t)> generateCombinations;
    generateCombinations = [&](size_t index)
    {
        // Base case: we have reached the end of the input elements,
        // so we add the current combination to the result
        if (index == elements.size())
        {
            combinations.push_back(currentCombination);
            return;
        }
        // Recursive case: we add the element at the current index to
        // the current combination and generate combinations starting
        // at the next index
        currentCombination.push_back(elements[index]);
        generateCombinations(index + 1);
        currentCombination.pop_back();
        // Recursive case: we skip the element at the current index and
        // generate combinations starting at the next index
        generateCombinations(index + 1);
    };
    // Start generating combinations at the first index
    generateCombinations(0);
    return combinations;
}

std::vector<cds::Residue*> residueCombinator::generateResidueCombinations(cds::Residue* starterResidue)
{
    // Find all positions that can be substituted
    std::vector<std::string> atomNumbers =
        residueCombinator::selectAllAtomsThatCanBeSubstituted(starterResidue->getAtoms());
    // Create the combinations of the numbers
    std::vector<std::vector<std::string>> numberCombinations = residueCombinator::getCombinations(atomNumbers);
    // Create a prep residue for each of the combinations, copying and modifying the original.
    std::vector<cds::Residue*> glycamResiduePermutations;
    for (auto& combination : numberCombinations)
    {
        glycamResiduePermutations.emplace_back(new cds::Residue(*starterResidue));
        cds::Residue* newResidue = glycamResiduePermutations.back();
        std::string delimiter    = "";
        std::stringstream numbersAsString;
        for (auto& atomNumber : combination)
        {
            numbersAsString << delimiter << atomNumber;
            delimiter           = ",";
            cds::Atom* hydrogen = newResidue->FindAtom("H" + atomNumber + "O");
            cds::Atom* oxygen   = newResidue->FindAtom("O" + atomNumber);
            oxygen->setCharge(oxygen->getCharge() + hydrogen->getCharge() - 0.194);
            newResidue->deleteAtom(hydrogen);
            oxygen->setType("Os");
        }
        std::string residueName = GlycamMetadata::GetGlycam06ResidueLinkageCode(numbersAsString.str());
        residueName             += starterResidue->getName().substr(1);
        newResidue->setName(residueName);
        std::cout << numbersAsString.str() << ": " << newResidue->getName() << "\n";
    }
    return glycamResiduePermutations;
}
