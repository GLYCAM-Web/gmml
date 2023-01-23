#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_

namespace cdsCondensedSequence
{
struct GraphVizDotConfig
{
    GraphVizDotConfig(std::string svg_directory_path = "/programs/gems/gmml/includes/MolecularMetadata/Sugars/SNFG_Symbol_Images/") :  svg_directory_path_(svg_directory_path) {}
    bool show_config_labels_ = true;
    bool show_edge_labels_ = false;
    bool show_position_labels_ = true;
    int dpi_ = 72;
    std::string file_name_ = "oligosaccharide.dot";
    std::string svg_directory_path_ = "/programs/gems/gmml/includes/MolecularMetadata/Sugars/SNFG_Symbol_Images/";
};
} // namespace
#endif /* INCLUDES_INPUTSET_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_ */
