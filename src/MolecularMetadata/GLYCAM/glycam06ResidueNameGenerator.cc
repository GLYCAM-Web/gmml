#include "includes/MolecularMetadata/GLYCAM/glycam06ResidueNameGenerator.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06LinkageCodes.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <locale> // for isLower()

namespace gmml
{
    namespace MolecularMetadata
    {
        namespace GLYCAM
        {
            std::string Glycam06ResidueNameGenerator(std::string linkages, std::string isomer, std::string inputResName,
                                                     std::string ringType, std::string residueModifier,
                                                     std::string configuration)
            {
                /* 	Example inputs:
                        linkages: "2,3" , "1" , "Terminal" , "4,7"
                        isomer: Always D or L
                        inputResName: "Glc" , "Neu" , "Ido"
                        ringType: Always f or p
                        residueModifier: "NAc" , "5Ac" , "A"
                        configuration: Always a or b

                        Example output:
                        UYB, 0GA etc. See glycam naming / nomenclature.
                */
                // Link code e.g. 0, 1, 2, W, Z etc
                std::string inputs = "\nInputs:\nlinkages: " + linkages + "\nisomer: " + isomer +
                                     "\ninputResName: " + inputResName + "\nringType: " + ringType +
                                     "\nresidueModifier: " + residueModifier + "\nconfiguration: " + configuration;
                // gmml::log(__LINE__, __FILE__, gmml::INF, inputs);

                std::string linkCode = ""; // Can be empty for e.g. ROH
                if (!linkages.empty())
                {
                    Glycam06LinkageCodesLookupContainer linkageCodeLookup;
                    linkCode = linkageCodeLookup.GetCodeForLinkages(linkages);
                    if (linkCode.empty())
                    {
                        std::string message = "No linkage code found in GMML metadata for a carbohydrate residue with "
                                              "other residues attached at these positions: " +
                                              linkages + "\nCheck these inputs for mistakes: " + inputs;
                        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
                        throw std::runtime_error(message);
                    }
                }

                // Configuration Code i.e. A/B/U/D
                std::string configurationCode = configuration; // I guess this is an ok default
                if ((ringType == "f") && (configuration == "a"))
                {
                    configurationCode = "D";
                }
                else if ((ringType == "f") && (configuration == "b"))
                {
                    configurationCode = "U";
                }
                else if ((ringType == "p") && (configuration == "a"))
                {
                    configurationCode = "A";
                }
                else if ((ringType == "p") && (configuration == "b"))
                {
                    configurationCode = "B";
                }

                // Residue Code e.g. G, U, A, KN, ROH, NLN
                Glycam06ResidueNamesToCodesLookupContainer ResidueCodeLookup;
                std::string residueCode = ResidueCodeLookup.GetCodeForResidue(inputResName + residueModifier);
                if (residueCode.empty())
                {
                    residueCode =
                        ResidueCodeLookup.GetCodeForResidue(inputResName + ringType + residueModifier + configuration);
                }
                if (residueCode.empty())
                {
                    residueCode = ResidueCodeLookup.GetCodeForResidue(isomer + inputResName + ringType +
                                                                      residueModifier + configuration);
                }
                if (residueCode.empty())
                {
                    residueCode =
                        ResidueCodeLookup.GetCodeForResidue(isomer + inputResName + ringType + residueModifier);
                }
                if (residueCode.empty())
                {
                    std::string message = "Cannot create 3D structure as no GLYCAM residue code was found in the GMML "
                                          "metadata for residue: " +
                                          isomer + inputResName + ringType + residueModifier + configuration;
                    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
                    throw std::runtime_error(message);
                }
                if (residueCode.size() > 1)
                {
                    configurationCode = ""; // Make it empty, it will implied by residueCode
                }
                // D vs L sugars. Residue code will be lowercase for L sugars
                if (isomer == "L") //&& (residueCode.size() == 1))
                {                  // Update: I think it should always lower, even for e.g. YS, KO, etc,
                    residueCode.at(0) = std::tolower(residueCode.at(0));
                }
                // ConfigurationCode may be empty.
                // gmml::log(__LINE__, __FILE__, gmml::INF,
                //          ("Returning: " + linkCode + residueCode + configurationCode + "\n"));
                return (linkCode + residueCode + configurationCode);
            }
        } // namespace GLYCAM
    }     // namespace MolecularMetadata
} // namespace gmml
