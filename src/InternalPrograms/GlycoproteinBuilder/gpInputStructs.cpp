#include "includes/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/InternalPrograms/functionsForGMML.hpp"
#include "includes/InternalPrograms/io.hpp"
#include "includes/CodeUtils/logging.hpp"

GlycoproteinBuilderInputs GPInputs::readGPInputFile(std::string workingDirectory, std::string inputFileName)
{
    std::ifstream infile (workingDirectory + inputFileName);
    if (!infile)
    {
        std::string message = "Uh oh, input file: " + workingDirectory + inputFileName + ", could not be opened for reading!\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw(message);
    }
    GlycoproteinBuilderInputs gpInputs;
    gpInputs.workingDirectory_ = workingDirectory; // Just to make sure. Having two like  directory// doesn't matter.
    while (infile) // While there's still stuff left to read
    {
        std::string strInput;
        getline(infile, strInput);
        if(gmml::startsWith(strInput, "Protein:"))
        {
        	gpInputs.substrateFileName_ = split(strInput, ':').at(1);
        }
        if(gmml::startsWith(strInput, "NumberOfOutputStructures:"))
        {
        	gpInputs.number3DStructures_ = split(strInput, ':').at(1);
        }
        if(gmml::startsWith(strInput, "prepFileLocation:"))
        {
        	gpInputs.prepFileLocation_ = split(strInput, ':').at(1);
        }
        if(gmml::startsWith(strInput, "maxThreads:"))
        {
        	gpInputs.maxThreads_ = split(strInput, ':').at(1);
        }
        if(gmml::startsWith(strInput, "persistCycles:"))
        {
        	gpInputs.persistCycles_ = split(strInput, ':').at(1);
        }
        if(gmml::startsWith(strInput, "overlapTolerance:"))
        {
        	gpInputs.overlapTolerance_ = split(strInput, ':').at(1);
        }
        if(gmml::startsWith(strInput, "isDeterministic:"))
        {
        	gpInputs.isDeterministic_ = split(strInput, ':').at(1);
        }
        if(strInput == "ProteinResidue, GlycanInputType, GlycanName:")
        {
            std::string tempBuffer; //  Temporarily holds whatever getline() finds on the line;
            getline(infile, tempBuffer);
            while(tempBuffer != "END")
            {
                std::vector<std::string> splitLine = split(tempBuffer, '|');
                gpInputs.glycositesInputVector_.emplace_back(splitLine.at(0), splitLine.at(1), splitLine.at(2));
                getline(infile, tempBuffer);
            }
        }
    }
    return gpInputs;
}



