#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <fstream>

using glycoprotein::GlycoproteinBuilderInputs;

GlycoproteinBuilderInputs glycoprotein::readGPInputFile(std::string inputFileName)
{
    //    std::cout << "About to read " << inputFileName << std::endl << std::flush;
    std::ifstream infile(inputFileName);
    if (!infile)
    {
        std::string message = "Uh oh, input file: " + inputFileName + ", could not be opened for reading!\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    GlycoproteinBuilderInputs gpInputs;
    while (infile) // While there's still stuff left to read
    {
        std::string strInput;
        getline(infile, strInput);
        if (codeUtils::startsWith(strInput, "Protein:"))
        {
            gpInputs.substrateFileName_ = codeUtils::split(strInput, ':').at(1);
        }
        if (codeUtils::startsWith(strInput, "NumberOfOutputStructures:"))
        {
            gpInputs.number3DStructures_ = std::stoi(codeUtils::split(strInput, ':').at(1));
        }
        if (codeUtils::startsWith(strInput, "maxThreads:"))
        {
            gpInputs.maxThreads_ = std::stoi(codeUtils::split(strInput, ':').at(1));
        }
        if (codeUtils::startsWith(strInput, "persistCycles:"))
        {
            gpInputs.persistCycles_ = std::stoi(codeUtils::split(strInput, ':').at(1));
        }
        if (codeUtils::startsWith(strInput, "overlapTolerance:"))
        {
            gpInputs.overlapTolerance_ = std::stoi(codeUtils::split(strInput, ':').at(1));
        }
        if (codeUtils::startsWith(strInput, "isDeterministic:"))
        { // variable = (condition) ? expressionTrue : expressionFalse;
            gpInputs.isDeterministic_ = (codeUtils::split(strInput, ':').at(1) == "true") ? true : false;
        }
        if (codeUtils::startsWith(strInput, "skipMDPrep:"))
        {
            gpInputs.skipMDPrep_ = (codeUtils::split(strInput, ':').at(1) == "true") ? true : false;
        }
        if (strInput == "ProteinResidue, GlycanName:")
        {
            std::string tempBuffer; //  Temporarily holds whatever getline() finds on the line;
            getline(infile, tempBuffer);
            while (tempBuffer != "END")
            {
                std::vector<std::string> splitLine = codeUtils::split(tempBuffer, '|');
                gpInputs.glycositesInputVector_.emplace_back(splitLine.at(0), splitLine.at(1));
                getline(infile, tempBuffer);
            }
        }
    }
    //    std::cout << "Reading input file complete, just making a quick check\n" << std::flush;
    if (gpInputs.glycositesInputVector_.empty())
    {
        throw std::runtime_error(
            "Error reading from gpInput file, no glycosites requested. Perhaps your formatting is incorrect.\n");
    }
    return gpInputs;
}
