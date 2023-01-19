#include "includes/CentralDataStructure/Readers/Prep/prepPermutator.hpp"

std::vector<std::string> prep::selectAllAtomsThatCanBeSubstituted(std::vector<cds::Atom*> queryAtoms)
{
    std::vector<std::string> foundNames;
    for(auto & atom : queryAtoms )
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

//ChatGPT made this:
std::vector<std::vector<std::string>> prep::getCombinations(const std::vector<std::string>& elements)
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


std::vector<prep::PrepResidue*> prep::generatePrepPermuations(PrepResidue* starterResidue)
{
    std::vector<PrepResidue*> prepPermutations;
    // Find all positions that can be substituted

    std::cout << "Finding hydroxyl oxygen numbers:" << std::endl;
    std::vector<std::string> atomNumbers = prep::selectAllAtomsThatCanBeSubstituted(starterResidue->getAtoms());
        // Create the combinations of the numbers
    std::cout << "Create combos:" << std::endl;
    std::vector<std::vector<std::string>> numberCombinations = prep::getCombinations(atomNumbers);
    // Create a prep residue for each of the combinations, copying and modifying the original.
    for(auto & combination : numberCombinations)
    {
        std::cout << "Creating a new residue" << std::endl;
        prepPermutations.emplace_back(new PrepResidue(*starterResidue));
        std::cout << "Fin" << std::endl;
        PrepResidue* newResidue = prepPermutations.back();
        std::cout << "Fin" << std::endl;
        for(auto & atomNumber : combination)
        {
            std::cout << "atomNumber: " << atomNumber  << std::endl;
            cds::Atom* hydrogen = newResidue->FindAtom("H" + atomNumber + "O");
            cds::Atom* oxygen = newResidue->FindAtom("O" + atomNumber);
            oxygen->setCharge(oxygen->getCharge() + hydrogen->getCharge() - 0.194);
            newResidue->deleteAtom(hydrogen);
            oxygen->setType("Os");
        }
        std::string residueName = ""; // // Use the class in #include "includes/MolecularMetadata/GLYCAM/glycam06LinkageCodes.hpp" to figure out the code
        newResidue->setName(residueName);
    }
    return prepPermutations;
}
