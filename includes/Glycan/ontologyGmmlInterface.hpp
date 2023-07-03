#ifndef GMML_INCLUDES_GLYCAN_ONTOLOGYGMMLINTERFACE_HPP
#define GMML_INCLUDES_GLYCAN_ONTOLOGYGMMLINTERFACE_HPP

#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include <sstream>

namespace Ontology
{
    void PrintOntology(std::stringstream& ont_stream, const pdb::PdbFile& pdbFile);
}

#endif // GMML_INCLUDES_GLYCAN_ONTOLOGYGMMLINTERFACE_HPP
