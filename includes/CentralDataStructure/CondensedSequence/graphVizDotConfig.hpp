#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_

#include <filesystem>

// unnamed namespace so we can grab the path we care about
namespace
{
    std::string getSNFGPath()
    {
        std::filesystem::path snfgSVGCollectionPath = __FILE__;
        snfgSVGCollectionPath                       = snfgSVGCollectionPath.parent_path().parent_path().parent_path();
        snfgSVGCollectionPath                       /= "MolecularMetadata/Sugars/SNFG_Symbol_Images/";
        return snfgSVGCollectionPath.string();
    }

    static const std::string snfgSymbolsPath = getSNFGPath();
} // namespace

namespace cdsCondensedSequence
{

    // To help intelligently set the path of all our svg dudes

    struct GraphVizDotConfig
    {
        GraphVizDotConfig(std::string svg_directory_path = snfgSymbolsPath) : svg_directory_path_(svg_directory_path)
        {}

        bool show_config_labels_   = true;
        bool show_edge_labels_     = false;
        bool show_position_labels_ = true;
        int dpi_                   = 72;
        std::string file_name_     = "oligosaccharide.dot";

        std::string svg_directory_path_ = snfgSymbolsPath;
    };

} // namespace cdsCondensedSequence
#endif /* INCLUDES_INPUTSET_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_ */
