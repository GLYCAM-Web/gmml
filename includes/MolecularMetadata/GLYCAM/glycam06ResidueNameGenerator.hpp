#ifndef GLYCAM06_RESIDUE_NAME_GENERATOR_HPP
#define GLYCAM06_RESIDUE_NAME_GENERATOR_HPP

#include <string>
#include <map>
#include <vector>

namespace gmml
{
namespace MolecularMetadata
{
namespace GLYCAM
{
std::string Glycam06ResidueNameGenerator(std::string linkages, char isomer, std::string inputResName, char ringType, std::string residueModifier, char configuration);
} // close namespace
} // close namespace
} // close namespace

#endif // GLYCAM06_RESIDUE_NAME_GENERATOR_HPP