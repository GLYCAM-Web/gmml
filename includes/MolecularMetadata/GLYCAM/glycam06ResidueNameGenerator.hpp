#ifndef GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUENAMEGENERATOR_HPP
#define GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUENAMEGENERATOR_HPP

#include <string>
#include <map>
#include <vector>

namespace gmml
{
namespace MolecularMetadata
{
namespace GLYCAM
{
std::string Glycam06ResidueNameGenerator(std::string linkages, std::string isomer, std::string inputResName, std::string ringType, std::string residueModifier, std::string configuration);
} // close namespace
} // close namespace
} // close namespace

#endif // GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUENAMEGENERATOR_HPP