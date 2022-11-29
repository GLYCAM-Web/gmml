#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_LIBRARYFILE_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_LIBRARYFILE_HPP_

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Readers/LibraryResidue.hpp"

#include <sstream>

namespace lib
{

class LibraryFile : public cds::Molecule
{
public:
    LibraryFile();
    LibraryFile(const std::string &filePath);
private:
    std::stringstream ExtractUnitSection(std::ifstream& inputFileStream, const std::string unitName);
    void ParseInFileStream(std::ifstream& inputFileStream);
};

} // namespace
#endif /* INCLUDES_CENTRALDATASTRUCTURE_READERS_LIBRARYFILE_HPP_ */
