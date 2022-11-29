#include "includes/CentralDataStructure/Readers/LibraryFile.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp" // RemoveSpaces/Quotes
#include <fstream> // std::ifstream

namespace lib
{
LibraryFile::LibraryFile(const std::string &filePath)
{
    codeUtils::ensureFileExists(filePath);
    std::ifstream fileStream(filePath);
    if(fileStream.fail())
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Could not open this file: " + filePath);
        throw std::runtime_error("PdbFile constructor could not open this file: " + filePath);
    }
    gmml::log(__LINE__, __FILE__,  gmml::INF, "File opened: " + filePath + ". Ready to parse!");
    this->ParseInFileStream(fileStream);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Finished parsing " + filePath);
}

// Private
void LibraryFile::ParseInFileStream(std::ifstream& inputFileStream)
{
    std::string line;
    // Skip any blank lines at the beginning of the file
    while(line.empty() || line.front() != '!')
    {
        getline(inputFileStream, line);
    }
    // Read in the array of residue names
    //std::cout << "First readable line is " << line << "\n";
    std::vector<std::string> residueNames;
    if(line.find("index") != std::string::npos)
    {
        while( getline(inputFileStream, line) && line.front() != '!')
        {
            codeUtils::RemoveQuotes(line);
            codeUtils::RemoveSpaces(line);
            residueNames.push_back(line);
            //
        }
    }
    if(line.find("atoms") == std::string::npos) // If first line passed in for a residue isn't the atoms table, Freak out.
    {
        std::string message = "Error reading library file, I expected the !entry line of an atoms table, but got this: " + line;
        throw std::runtime_error(message);
    }
    // Iterate on residue names
    for(auto & residueName : residueNames)
    {   // Process the atom section of the file for the corresponding residue

        //std::cout << "Starting to read residue " << residueName << " with line:\n" << line << "\n";
        std::stringstream residueStream;
        residueStream = this->ExtractUnitSection(inputFileStream, residueName);
        //std::cout << "Residue stream is:\n" << residueStream.str() << "\n fin. " << std::endl;
        this->addResidue(std::make_unique<LibraryResidue>(residueStream, residueName));
    }
}

std::stringstream LibraryFile::ExtractUnitSection(std::ifstream& inputFileStream, const std::string unitName)
{
    std::stringstream extractedSection;
    std::string entryLineStart = "!entry." + unitName + ".unit";
    std::string line;
    while ((std::getline(inputFileStream, line)))
    {
        if ( line.find("!entry") != std::string::npos && line.find(entryLineStart) == std::string::npos )
        { // If line is an "entry" line that doesn't match this unit, time to leave.
            //std::cout << "I have never seen this line before in my life: " << line << std::endl;
            return extractedSection;
        }
        extractedSection << line << std::endl;
    }
    return extractedSection;
}

} // namespace
