#include "includes/CentralDataStructure/Editors/glycamResidueCombinator.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" //calculateCoordinateFromInternalCoords

void residueCombinator::removeHydroxyHydrogen(cds::Residue& queryResidue, const std::string hydrogenNumber)
{
    cds::Atom* hydrogen = queryResidue.FindAtom("H" + hydrogenNumber + "O");
    cds::Atom* oxygen   = queryResidue.FindAtom("O" + hydrogenNumber);
    if (hydrogen == nullptr || oxygen == nullptr)
    {
        std::string message =
            "Cannot find appropriately named atoms in residue. Glycam combinations cannot be created. Oxygen should be "
            "named e.g. O2 and not 2O. Hydrogen to be substituted should be H2O and not HO2. Both must be present. "
            "This may turn into a fatal issue for the atom numbered: " +
            hydrogenNumber + " in residue: " + queryResidue.getName();
        gmml::log(__LINE__, __FILE__, gmml::WAR, message);
        return;
    }
    oxygen->setCharge(oxygen->getCharge() + hydrogen->getCharge() - 0.194);
    queryResidue.deleteAtom(hydrogen);
    oxygen->setType("Os");
}

std::vector<std::string> residueCombinator::selectAllAtomsThatCanBeSubstituted(const cds::Residue& queryResidue)
{
    std::vector<std::string> foundNames;
    std::string delimiter = "";
    for (auto& atom : queryResidue.getAtoms())
    { // if a hydroxyl with a digit in the second position of atom name. e.g. O2, not OHG, and not O1A.
        if (atom->getType() == "Oh" && atom->getNumberFromName() != 0 &&
            isdigit(atom->getName().back())) // If a hydroxyl like O2
        {
            foundNames.push_back(std::to_string(atom->getNumberFromName()));
            std::cout << delimiter << atom->getName();
            delimiter = ",";
        }
    }
    std::cout << "\n";
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

void generateResidueCombination(std::vector<cds::Residue*>& glycamResidueCombinations,
                                const std::vector<std::string> numberCombination, const cds::Residue& templateResidue)
{
    cds::Residue* newResidue = glycamResidueCombinations.emplace_back(new cds::Residue(templateResidue));
    // cds::Residue* newResidue = glycamResidueCombinations.back();
    std::string delimiter    = "";
    std::stringstream numbersAsString;
    for (auto& atomNumber : numberCombination)
    {
        numbersAsString << delimiter << atomNumber;
        delimiter = ",";
        residueCombinator::removeHydroxyHydrogen(*newResidue, atomNumber);
        // Set the tail. All can go to the end of the residue. Don't think order matters.
        Atom* atom = newResidue->FindAtom("O" + atomNumber);
        newResidue->moveAtomToLastPosition(atom);
    }
    std::string residueName = GlycamMetadata::GetGlycam06ResidueLinkageCode(numbersAsString.str());
    if (residueName.empty())
    { // Now we have the need for a new residue nomenclature to kick in.
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "No linkage code found for possible combo: " + numbersAsString.str() + " in residue " +
                      templateResidue.getName());
        glycamResidueCombinations.pop_back(); // This should be rare/never?
    }
    else
    {
        residueName += templateResidue.getName().substr(1);
        newResidue->setName(residueName);
        std::cout << numbersAsString.str() << ": " << newResidue->getName() << "\n";
    }
}

// Problem: The input residue will either contain an anomeric oxygen or not.
// Must generate the version that's missing as both are required.
// Combinations that include the anomeric position should include an anomeric oxygen (eg 1GA).
// Combinations that do not include the anomeric position should not (eg 0GA).

void residueCombinator::generateResidueCombinations(std::vector<cds::Residue*>& glycamResidueCombinations,
                                                    const cds::Residue* starterResidue)
{
    // First generate both versions of the residue; with and without anomeric oxygen.
    cds::Residue residueWithoutAnomericOxygen = *starterResidue;
    cds::Residue residueWithAnomericOxygen    = *starterResidue;
    Atom* anomer             = cdsSelections::guessAnomericAtomByInternalNeighbors(starterResidue->getAtoms());
    std::string anomerNumber = std::to_string(anomer->getNumberFromName());
    std::vector<Atom*> anomerNeighbors = anomer->getNeighbors();
    auto isTypeHydroxy                 = [](Atom*& a) // Lamda function for the  std::find;
    {
        return (a->getType() == "Oh");
    };
    const auto anomericOxygen = std::find_if(anomerNeighbors.begin(), anomerNeighbors.end(), isTypeHydroxy);
    if (anomericOxygen != anomerNeighbors.end())
    {
        std::cout << "Anomeric Oxygen Found\n";
        residueCombinator::removeHydroxyHydrogen(residueWithAnomericOxygen, anomerNumber);
        residueWithoutAnomericOxygen = residueWithAnomericOxygen;
        residueWithoutAnomericOxygen.deleteAtom(*anomericOxygen);
        // ToDo CHARGE?
    }
    else // Ok then grow the anomeric oxygen for the residueWithAnomericOxygen
    {
        std::cout << "No Anomeric Oxygen Found in Template\n";
        std::vector<Coordinate*> threeNeighbors;
        for (auto& neighbor : anomer->getNeighbors())
        {
            threeNeighbors.push_back(neighbor->getCoordinate());
        }
        Coordinate newOxygenCoordinate =
            cds::CreateCoordinateForCenterAwayFromNeighbors(anomer->getCoordinate(), threeNeighbors, 1.4);
        Atom* newAnomericOxygen = residueWithAnomericOxygen.addAtomToFront(
            std::make_unique<cds::Atom>("O" + anomerNumber, newOxygenCoordinate));
        newAnomericOxygen->setCharge(-0.388);
        newAnomericOxygen->setType("Os");
    }
    // Find all positions that can be substituted, ignore the anomer.
    std::vector<std::string> atomNumbers =
        residueCombinator::selectAllAtomsThatCanBeSubstituted(residueWithoutAnomericOxygen);
    // Create the combinations of the numbers
    std::vector<std::vector<std::string>> numberCombinations = residueCombinator::getCombinations(atomNumbers);
    // Create a residue for each of the combinations, copying and modifying the original.
    for (auto& combination : numberCombinations)
    {
        generateResidueCombination(glycamResidueCombinations, combination, residueWithoutAnomericOxygen);
        // ToDo Activate the below when we can handle combinations with anomeric positions.
        // combination.push_back(anomerNumber);
        // generateResidueCombination(glycamResidueCombinations, combination, residueWithAnomericOxygen);
    }
    // Handle the anomeric position separately. No combinations allowed with this position.
    cds::Residue* newResidue = glycamResidueCombinations.emplace_back(new cds::Residue(residueWithAnomericOxygen));
    std::string residueName  = anomerNumber;
    residueName              += starterResidue->getName().substr(1);
    newResidue->setName(residueName);
    std::cout << "Added " << residueName << "\n";
    // Write out the 0.. version:
    newResidue  = glycamResidueCombinations.emplace_back(new cds::Residue(residueWithoutAnomericOxygen));
    residueName = "0";
    residueName += starterResidue->getName().substr(1);
    newResidue->setName(residueName);
    std::cout << "Added " << residueName << "\n";
    return;
}
