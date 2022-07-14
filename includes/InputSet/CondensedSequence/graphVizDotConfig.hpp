#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_

#include <string>

namespace CondensedSequence
{
struct GraphVizDotConfig
{
    bool show_config_labels_;
    bool show_edge_labels_;
    bool show_position_labels_;
    int dpi_;
    std::string file_name_;
    std::string svg_directory_path_;
    GraphVizDotConfig()
    {
        this->show_edge_labels_ = false;
        this->show_config_labels_ = true;
        this->show_position_labels_ = true;
        this->dpi_ = 72;
        this->svg_directory_path_ = "/programs/gems/gmml/includes/MolecularMetadata/Sugars/SNFG_Symbol_Images/";
        this->file_name_ = "oligosaccharide.dot";
    }
};
} // namespace
#endif /* INCLUDES_INPUTSET_CONDENSEDSEQUENCE_GRAPHVIZDOTCONFIG_HPP_ */
