#ifndef GLYCAM06_RESIDUE_NAMES_TYPES_META_HPP
#define GLYCAM06_RESIDUE_NAMES_TYPES_META_HPP

/** \file:  includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp
 * GLYCAM06 metadata for residues
 *
 * This file was generated automatically on:
 *     Thu Sep 13 11:01:29 EDT 2018
 *
 * by a script named:
 *     1_make_glycam06_residue_tag_lookup.bash
 *
 * The script was begun on 16 June 2018 by BLFoley and
 * can be found in:
 *     dat/MolecularMetadata/scripts/GLYCAM/GLYCAM06/ResidueTags
 *
 * See that and associated scripts for more information.
 */

#include <string>
#include <map>
#include <vector>

namespace gmml
{
namespace MolecularMetadata
{
namespace GLYCAM
{
/**
         *   Glycam06NamesToTypesLookupMap
         *
         *   The first string is the name-code for a residue.  It is typically
         *   three or four characters long.
         *
         *     Examples:  OME, 0GA, WYB, ZOLT YuAP
         *
         *   The second string is a type or attribute appropriate to that residue.
         *
         *     Examples:  monosaccharide, derivative, aglycon, protonated, pyranose
         *
         *
         */

class Glycam06NamesToTypesLookupContainer
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
    * Default constructor
    */
    Glycam06NamesToTypesLookupContainer(); // Calls an initializer?

    //////////////////////////////////////////////////////////
    //                         TYPEDEFS                     //
    //////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    //                      QUERY FUNCTIONS                 //
    //////////////////////////////////////////////////////////

    inline std::vector<std::string> GetTypesForResidue(std::string query)
    {
        std::vector<std::string> matching_types;
        // Iterate over the multimap using range based for loop
        for (std::pair<std::string, std::string> elem : glycam06NamesToTypesLookupMap_)
        {
            if (elem.first.compare(query)==0)
            {
                matching_types.push_back(elem.second);
            }
            //std::cout << elem.first << " :: " << elem.second << std::endl;
        }
        return matching_types;
    }
private:
    std::multimap<std::string, std::string> glycam06NamesToTypesLookupMap_;
};
} // close namespace
} // close namespace
} // close namespace

#endif // GLYCAM06_RESIDUE_NAMES_TYPES_META_HPP
