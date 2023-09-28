#include "includes/MolecularMetadata/proteinBonding.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <unordered_map>
#include <stdexcept>

const std::vector<std::pair<std::string, std::string>>& biology::getBackboneBonding()
{
    static const std::vector<std::pair<std::string, std::string>> backboneBonding = {
        { "N", "CA"},
        {"CA",  "C"},
        { "C",  "O"}
    };
    return backboneBonding;
}

const std::string biology::equivalentResidueNameLookup(const std::string queryResidueName)
{
    static const std::unordered_map<std::string, std::string> proteinNameMap = {
        {"HIE", "HIS"},
        {"HID", "HIS"},
        {"HIP", "HIS"},
        {"CYX", "CYS"},
        {"CYM", "CYS"},
        {"NLN", "ASN"},
        {"OLS", "SER"},
        {"OLT", "THR"},
        {"OLY", "TYR"},
        {"ARN", "ARG"},
        {"ASH", "ASP"},
        {"GLH", "GLU"},
        {"LYN", "LYS"},
    };
    if (proteinNameMap.find(queryResidueName) == proteinNameMap.end())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Found nothing in weirdo lookup map for " + queryResidueName);
        return queryResidueName;
    }
    gmml::log(__LINE__, __FILE__, gmml::INF,
              "Found a conversion for " + queryResidueName + " : " + proteinNameMap.at(queryResidueName));
    return proteinNameMap.at(queryResidueName);
}

const std::vector<std::pair<std::string, std::string>>& biology::getSidechainBonding(std::string queryResidueName)
{
    static const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> sidechainBonding = {
        {"ALA",                                                                        {{"CA", "CB"}}               },
        {"ARG", {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "NE"}, {"NE", "CZ"}, {"CZ", "NH1"}, {"CZ", "NH2"}}},
        {"ASP",                                           {{"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "OD2"}}},
        {"ASN",                                           {{"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "ND2"}}},
        {"CYS",                                                                         {{"CA", "CB"}, {"CB", "SG"}}},
        {"GLU",                             {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "OE2"}}},
        {"GLN",                             {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "NE2"}}},
        {"GLY",                                                                                                   {}}, // Can we handle this nicely?
        {"HIS",
         {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD2"}, {"CD2", "NE2"}, {"CG", "ND1"}, {"ND1", "CE1"}, {"CE1", "NE2"}} },
        {"ILE",                                         {{"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}, {"CG1", "CD1"}}},
        {"LEU",                                           {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}}},
        {"LYS",                               {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "CE"}, {"CE", "NZ"}}},
        {"MET",                                             {{"CA", "CB"}, {"CB", "CG"}, {"CG", "SD"}, {"SD", "CE"}}},
        {"PHE",
         {{"CA", "CB"},
         {"CB", "CG"},
         {"CG", "CD1"},
         {"CD1", "CE1"},
         {"CE1", "CZ"},
         {"CG", "CD2"},
         {"CD2", "CE2"},
         {"CE2", "CZ"}}                                                                                             },
        {"PRO",                                              {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "N"}}},
        {"SER",                                                                         {{"CA", "CB"}, {"CB", "OG"}}},
        {"THR",                                                         {{"CA", "CB"}, {"CB", "OG1"}, {"CB", "CG2"}}},
        {"TRP",
         {{"CA", "CB"},
         {"CB", "CG"},
         {"CG", "CD1"},
         {"CG", "CD2"},
         {"CD1", "NE1"},
         {"CD2", "CE2"},
         {"CD2", "CE3"},
         {"NE1", "CE2"},
         {"CE2", "CZ2"},
         {"CE3", "CZ3"},
         {"CZ2", "CH2"},
         {"CZ3", "CH2"}}                                                                                            },
        {"TYR",
         {{"CA", "CB"},
         {"CB", "CG"},
         {"CG", "CD1"},
         {"CD1", "CE1"},
         {"CE1", "CZ"},
         {"CG", "CD2"},
         {"CD2", "CE2"},
         {"CE2", "CZ"},
         {"CZ", "OH"}}                                                                                              },
        {"VAL",                                                         {{"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}}},
        {"SEC",                                                                         {{"CA", "CB"}, {"CB", "SE"}}},
        {"MSE",                                             {{"CA", "CB"}, {"CB", "CG"}, {"CG", "SE"}, {"SE", "CE"}}},
    };
    queryResidueName =
        biology::equivalentResidueNameLookup(queryResidueName); // returns same string if didn't find conversion.
    if (sidechainBonding.find(queryResidueName) == sidechainBonding.end())
    {
        throw std::runtime_error("Oliver! Each query residue residue should be in biology::proteinResidueNames (you "
                                 "need to confirm this before calling this code), or if it's there then it needs an "
                                 "entry in this table. Git gud son. This was the problem residue: " +
                                 queryResidueName);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Found sidechain bonding info for " + queryResidueName);
    return (sidechainBonding.at(queryResidueName));
}
