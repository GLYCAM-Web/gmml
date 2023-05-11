#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/inputs.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <fstream>
#include <iostream>

using gmmlPrograms::WiggleToSiteInputs;

WiggleToSiteInputs::WiggleToSiteInputs(std::string inputFileName)
{
//    std::cout << "About to read " << inputFileName << std::endl << std::flush;
    std::ifstream infile (inputFileName);
    if (!infile)
    {
        std::string message = "Uh oh, input file: " + inputFileName + ", could not be opened for reading!\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    while (infile) // While there's still stuff left to read
    {
        std::string strInput;
        getline(infile, strInput);
        if(codeUtils::startsWith(strInput, "Substrate:"))
        {
            substrateFile_ = codeUtils::split(strInput, ':').at(1);
        }
        if(codeUtils::startsWith(strInput, "Carbohydrate:"))
        {
            carbohydrateSequence_ = codeUtils::split(strInput, ':').at(1);
        }
        if(codeUtils::startsWith(strInput, "SuperimpositionTargetResidue:"))
        {
            std::string inputPortion = codeUtils::split(strInput, ':').at(1);
            superimpositionTargetResidue_ = pdb::ResidueId(codeUtils::split(inputPortion, '_'));
//            std::cout << "superimpositionTargetResidue_ in input file is: " << superimpositionTargetResidue_ << std::endl;
        }
        if(codeUtils::startsWith(strInput, "WiggleTargetResidue:"))
        {
            std::string inputPortion = codeUtils::split(strInput, ':').at(1);
            wigglingTargetResidue_ = pdb::ResidueId(codeUtils::split(inputPortion, '_'));
//            std::cout << "wigglingTargetResidue_ in input file is: " << wigglingTargetResidue_ << std::endl;
        }
        if(codeUtils::startsWith(strInput, "TargetModelNumber:"))
        {
            substrateModelNumber_ = std::stoi(codeUtils::split(strInput, ':').at(1));
        }
        if(codeUtils::startsWith(strInput, "CarbohydrateSuperimpositionResidue:"))
        {
            carbohydrateSuperimpositionResidue_ = std::stoi(codeUtils::split(strInput, ':').at(1));
        }
        if(codeUtils::startsWith(strInput, "CarbohydrateWigglingResidue:"))
        {
            carbohydrateWigglingResidue_ = std::stoi(codeUtils::split(strInput, ':').at(1));
        }
        if(codeUtils::startsWith(strInput, "persistCycles:"))
        {
            persistCycles_ = std::stoi(codeUtils::split(strInput, ':').at(1));
        }
        if(codeUtils::startsWith(strInput, "isDeterministic:"))
        {
            if (codeUtils::split(strInput, ':').at(1) == "true")
            {
                isDeterministic_ = true;
            }
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Reading input file complete." );
}

std::string WiggleToSiteInputs::Print()
{
    std::stringstream ss;
    ss << "substrateFile_ : " << substrateFile_ << "\n";
    ss << "carbohydrateFile_ : " << carbohydrateSequence_ << "\n";
    ss << "superimpositionTargetResidue_ : " << superimpositionTargetResidue_ << "\n";
    ss << "carbohydrateSuperimpositionResidue_ : " << carbohydrateSuperimpositionResidue_ << "\n";
    ss << "wigglingTargetResidue_ : " << wigglingTargetResidue_ << "\n";
    ss << "carbohydrateWigglingResidue_ : " << carbohydrateWigglingResidue_ << "\n";
    ss << "persistCycles_ : " << persistCycles_ << "\n";
    ss << "isDeterministic_ : " << isDeterministic_ << "\n";
    return ss.str();
}



