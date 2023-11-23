#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/Editors/glycamResidueCombinator.hpp"
#include <fstream>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputPrepFile\n";
        std::cout << "Exmpl: " << argv[0] << "../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep";
        std::exit(1);
    }
    std::string inputFile = argv[1];
    std::cout << "Test 027 Input file is " << inputFile << "\n";
    // std::vector<std::string> residuesToLoadFromPrep = {"0GA", "0YB"};
    std::vector<std::string> residuesToLoadFromPrep = {
        {"0AA"}, {"0AB"}, {"0AD"}, {"0AU"}, {"0BA"}, {"0BB"}, {"0BC"}, {"0BD"}, {"0BU"}, {"0CA"}, {"0CB"}, {"0CD"},
        {"0CU"}, {"0DA"}, {"0DB"}, {"0DD"}, {"0DU"}, {"0EA"}, {"0EB"}, {"0FA"}, {"0FB"}, {"0GA"}, {"0GB"}, {"0GL"},
        {"0HA"}, {"0HB"}, {"0JA"}, {"0JB"}, {"0JD"}, {"0JU"}, {"0KA"}, {"0KB"}, {"0LA"}, {"0LB"}, {"0MA"}, {"0MB"},
        {"0NA"}, {"0NB"}, {"0OA"}, {"0OB"}, {"0PA"}, {"0PB"}, {"0PD"}, {"0PU"}, {"0QA"}, {"0QB"}, {"0RA"}, {"0RB"},
        {"0RD"}, {"0RU"}, {"0SA"}, {"0SB"}, {"0TA"}, {"0TB"}, {"0TV"}, {"0Tv"}, {"0UA"}, {"0UB"}, {"0VA"}, {"0VB"},
        {"0WA"}, {"0WB"}, {"0XA"}, {"0XB"}, {"0XD"}, {"0XU"}, {"0YA"}, {"0YB"}, {"0ZA"}, {"0ZB"}, {"0aA"}, {"0aB"},
        {"0aD"}, {"0aU"}, {"0bA"}, {"0bB"}, {"0bC"}, {"0bD"}, {"0bU"}, {"0cA"}, {"0cB"}, {"0cD"}, {"0cU"}, {"0dA"},
        {"0dB"}, {"0dD"}, {"0dU"}, {"0eA"}, {"0eB"}, {"0fA"}, {"0fB"}, {"0gA"}, {"0gB"}, {"0gL"}, {"0hA"}, {"0hB"},
        {"0jA"}, {"0jB"}, {"0jD"}, {"0jU"}, {"0kA"}, {"0kB"}, {"0lA"}, {"0lB"}, {"0mA"}, {"0mB"}, {"0nA"}, {"0nB"},
        {"0oA"}, {"0oB"}, {"0pA"}, {"0pB"}, {"0pD"}, {"0pU"}, {"0qA"}, {"0qB"}, {"0rA"}, {"0rB"}, {"0rD"}, {"0rU"},
        {"0sA"}, {"0sB"}, {"0tA"}, {"0tB"}, {"0tV"}, {"0tv"}, {"0uA"}, {"0uB"}, {"0vA"}, {"0vB"}, {"0wA"}, {"0wB"},
        {"0xA"}, {"0xB"}, {"0xD"}, {"0xU"}, {"0yA"}, {"0yB"}, {"0zA"}, {"0zB"}, {"0dR"}, {"045"}, {"0Yn"}, {"0YS"},
        {"0Ys"}, {"0yS"}, {"0ys"}, {"0Kn"}, {"0Ko"}, {"0KN"}, {"0KO"}};
    prep::PrepFile glycamPrepFile(inputFile, residuesToLoadFromPrep);
    std::vector<cds::Residue*> allTheResidues;
    allTheResidues.reserve(residuesToLoadFromPrep.size() * 10);
    for (auto& residue : glycamPrepFile.getResidues())
    {
        std::cout << "Generating those combos boi" << std::endl;
        std::vector<cds::Residue*> newResidues = residueCombinator::generateResidueCombinations(residue);
        newResidues.push_back(residue);
        std::ofstream outFileStream;
        std::string fileName = residue->getName() + ".lib";
        outFileStream.open(fileName.c_str());
        cds::WriteResiduesToOffFile(newResidues, outFileStream);
        outFileStream.close();
        allTheResidues.insert(allTheResidues.end(), newResidues.begin(), newResidues.end());
    }
    std::ofstream outFileStream;
    std::string fileName = "GLYCAM06k.lib";
    outFileStream.open(fileName.c_str());
    cds::WriteResiduesToOffFile(allTheResidues, outFileStream);
    outFileStream.close();
}
