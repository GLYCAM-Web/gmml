#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP_
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp" // serializeNumbers
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip> //std::precision

namespace cds
{
std::string getOffType(const cds::ResidueType queryType);
void WriteOffFileUnit(std::vector<cds::Residue*> residues, std::ostream& stream, const std::string unitName);
void WriteResiduesToOffFile(std::vector<cds::Residue*> residues, std::ostream& stream);
void WriteMoleculeToOffFile(cds::Molecule* molecule, std::ostream& stream, const std::string unitName);
void WriteAssemblyToOffFile(cds::Assembly* assembly, std::ostream& stream, const std::string unitName);
} // namespace
#endif
