#include "includes/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/InternalPrograms/functionsForGMML.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/strings.hpp"

GlycoproteinBuilderInputs GPInputs::readGPInputFile(std::string workingDirectory, std::string inputFileName)
{
    if (workingDirectory == "Default")
    {
        workingDirectory = codeUtils::Find_Program_workingDirectory();
    }
    std::ifstream infile (workingDirectory + inputFileName);
    if (!infile)
    {
        std::string message = "Uh oh, input file: " + workingDirectory + inputFileName + ", could not be opened for reading!\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    GlycoproteinBuilderInputs gpInputs;
    gpInputs.workingDirectory_ = workingDirectory; // Just to make sure. Having two like  directory// doesn't matter.
    while (infile) // While there's still stuff left to read
    {
        std::string strInput;
        getline(infile, strInput);
        if(codeUtils::startsWith(strInput, "Protein:"))
        {
        	gpInputs.substrateFileName_ = codeUtils::split(strInput, ':').at(1);
        }
        if(codeUtils::startsWith(strInput, "NumberOfOutputStructures:"))
        {
        	gpInputs.number3DStructures_ = codeUtils::split(strInput, ':').at(1);
        }
        if(codeUtils::startsWith(strInput, "prepFileLocation:"))
        {
        	gpInputs.prepFileLocation_ = codeUtils::split(strInput, ':').at(1);
        }
        if(codeUtils::startsWith(strInput, "maxThreads:"))
        {
        	gpInputs.maxThreads_ = codeUtils::split(strInput, ':').at(1);
        }
        if(codeUtils::startsWith(strInput, "persistCycles:"))
        {
        	gpInputs.persistCycles_ = codeUtils::split(strInput, ':').at(1);
        }
        if(codeUtils::startsWith(strInput, "overlapTolerance:"))
        {
        	gpInputs.overlapTolerance_ = codeUtils::split(strInput, ':').at(1);
        }
        if(codeUtils::startsWith(strInput, "isDeterministic:"))
        {
        	gpInputs.isDeterministic_ = codeUtils::split(strInput, ':').at(1);
        }
        if(strInput == "ProteinResidue, GlycanInputType, GlycanName:")
        {
            std::string tempBuffer; //  Temporarily holds whatever getline() finds on the line;
            getline(infile, tempBuffer);
            while(tempBuffer != "END")
            {
                std::vector<std::string> splitLine = codeUtils::split(tempBuffer, '|');
                gpInputs.glycositesInputVector_.emplace_back(splitLine.at(0), splitLine.at(1), splitLine.at(2));
                getline(infile, tempBuffer);
            }
        }
    }
    return gpInputs;
}



