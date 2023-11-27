#include "includes/CentralDataStructure/Editors/glycamResidueCombinator.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06LinkageCodes.hpp"
#include "includes/CodeUtils/logging.hpp"

std::vector<std::string> residueCombinator::selectAllAtomsThatCanBeSubstituted(std::vector<cds::Atom*> queryAtoms)
{
    std::vector<std::string> foundNames;
    std::string delimiter = "";
    for (auto& atom : queryAtoms)
    { // if a hydroxyl with a digit in the second position of atom name. e.g. O2, not OHG.
        if (atom->getType() == "Oh" && isdigit(atom->getName().at(1))) // If a hydroxyl.
        {
            foundNames.push_back(atom->getName().substr(1, 1));
            //            std::cout << delimiter << atom->getName();
            //            delimiter = ",";
        }
    }
    //    std::cout << "\nDone selecting atoms" << std::endl;
    return foundNames;
}

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
    combinations.pop_back(); // I don't want the last empty one.
    return combinations;
}

std::vector<cds::Residue*> residueCombinator::generateResidueCombinations(cds::Residue* starterResidue)
{
    std::vector<cds::Residue*> glycamResiduePermutations;
    // Find all positions that can be substituted
    std::vector<std::string> atomNumbers =
        residueCombinator::selectAllAtomsThatCanBeSubstituted(starterResidue->getAtoms());
    // Handle the special cases of reducing O1.
    auto found = std::find(atomNumbers.begin(), atomNumbers.end(), "1");
    if (found != std::end(atomNumbers))
    {
        atomNumbers.erase(found); // don't make combos with 1
        glycamResiduePermutations.emplace_back(new cds::Residue(*starterResidue));
        cds::Residue* newResidue = glycamResiduePermutations.back();
        newResidue->RemoveHydroxyHydrogen("1");
        newResidue->setName("1" + starterResidue->getName().substr(1));
    }
    // Create the combinations of the numbers
    std::vector<std::vector<std::string>> numberCombinations = residueCombinator::getCombinations(atomNumbers);
    // Create a prep residue for each of the combinations, copying and modifying the original.
    for (auto& combination : numberCombinations)
    {
        glycamResiduePermutations.emplace_back(new cds::Residue(*starterResidue));
        cds::Residue* newResidue = glycamResiduePermutations.back();
        std::string delimiter    = "";
        std::stringstream numbersAsString;
        for (auto& atomNumber : combination)
        {
            numbersAsString << delimiter << atomNumber;
            delimiter = ",";
            newResidue->RemoveHydroxyHydrogen(atomNumber);
        }
        std::string residueName = GlycamMetadata::GetGlycam06ResidueLinkageCode(numbersAsString.str());
        if (residueName.empty())
        { // Now we have the need for a new residue nomenclature to kick in.
            gmml::log(__LINE__, __FILE__, gmml::WAR,
                      "No linkage code found for possible combo: " + numbersAsString.str() + " in residue " +
                          starterResidue->getName());
        }
        else
        {
            residueName += starterResidue->getName().substr(1);
            newResidue->setName(residueName);
            // std::cout << numbersAsString.str() << ": " << newResidue->getName() << "\n";
        }
    }
    return glycamResiduePermutations;
}
