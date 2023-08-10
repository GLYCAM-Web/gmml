#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_

#include <filesystem>

// unnamed namespace so we can grab the path we care about
namespace
{
    // TODO: Do this not disgusting
    static const std::string snfgSymbolsPath =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path().string() +
        "/MolecularMetadata/Sugars/SNFG_Symbol_Images/";
} // namespace

namespace cdsCondensedSequence
{
    struct GraphVizDotConfig
    {
        GraphVizDotConfig(std::string svg_directory_path = snfgSymbolsPath) : svg_directory_path_(svg_directory_path)
        {}

        bool show_config_labels_        = true;
        bool show_edge_labels_          = false;
        bool show_position_labels_      = true;
        int dpi_                        = 72;
        std::string file_name_          = "oligosaccharide.dot";
        std::string svg_directory_path_ = snfgSymbolsPath;
    };

} // namespace cdsCondensedSequence
#endif /* INCLUDES_INPUTSET_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_ */
