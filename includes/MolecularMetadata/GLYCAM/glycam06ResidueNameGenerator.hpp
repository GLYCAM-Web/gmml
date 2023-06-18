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
            std::string Glycam06ResidueNameGenerator(std::string linkages, std::string isomer, std::string inputResName,
                                                     std::string ringType, std::string residueModifier,
                                                     std::string configuration);
        } // namespace GLYCAM
    }     // namespace MolecularMetadata
} // namespace gmml

#endif // GLYCAM06_RESIDUE_NAME_GENERATOR_HPP