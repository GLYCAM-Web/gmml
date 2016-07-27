#ifndef OLIGOSACCHARIDE_HPP
#define OLIGOSACCHARIDE_HPP

#include <vector>
#include <string>
#include <sstream>

#include "../MolecularModeling/assembly.hpp"
#include "../utils.hpp"

namespace Glycan
{
    struct Oligosaccharide
    {
            Monosaccharide* root_;                              /*!< The root sacchride of the oligosacchride >*/
            std::vector<Oligosaccharide*> child_oligos_;        /*!< The oligosacchrides that are attached to the current oligosaccharide >*/
            std::vector<std::string> child_oligos_linkages_;    /*!< The linkages between the current oligosaccharide and oligosacchrides that are attached to it. i.e. C4_13_4GB_?_2_?_?_1-O4_23_4GB_?_2_?_?_1-C1_24_0MA_?_3_?_?_1 >*/
            std::string oligosaccharide_name_;                  /*!< The complete name sequence of the oligosacchride >*/
            std::string oligosaccharide_linkages_;              /*!< The complete sequence of the linkages of the oligosacchride >*/
            std::string oligosaccharide_residue_linkages_;      /*!< The complete sequence of the residue linkages of the oligosacchride >*/
            std::string terminal_;                              /*!< The terminal residue name of the oligosacchride >*/

            /*! \fn
              * A function to print out the oligosacchride name and the linkages between its monosacchrides
              * Print out the information in a defined structure              
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout)
            {
                oligosaccharide_name_ = "";
                oligosaccharide_linkages_ = "";
                oligosaccharide_residue_linkages_ = "";
                int linkage_index = 0;
                std::stringstream ss;
                if(terminal_.compare("") != 0)
                {
                    if(root_->derivatives_map_.find("-1") == root_->derivatives_map_.end())
                        ss << "1-" << terminal_;
                    else
                        ss << "2-" << terminal_;
                }
                oligosaccharide_name_ = ss.str();
                bool is_cycle = false;
                this->GenerateNameLinkage(oligosaccharide_name_, oligosaccharide_linkages_, linkage_index, root_->mono_id, is_cycle);

                if(is_cycle)
                {
                    int end_index = oligosaccharide_name_.find_last_of('-');
                    std::vector<std::string> tokens = gmml::Split(oligosaccharide_name_, "-");

                    std::string sub_name = oligosaccharide_name_.substr(0, end_index);
                    std::stringstream new_name;
                    new_name << "[" << tokens.at(tokens.size() - 1).at(0) << sub_name << "-]";
                    oligosaccharide_name_ = new_name.str();
                }

                // Example oligo sequence LRhapa1-3LRhapa1-3DGlcpNAcb1-OME
                // Example oligo linkages
                //{2}RAM(401_A)C1-RAE(402_A)C3, Glycosidic linkage: RAE(402_A)O3
                //{1}RAE(402_A)C1-MAG(403_A)C3, Glycosidic linkage: MAG(403_A)O3
                //For each line of linkages the info of the first residue (e.g. RAE(402_A)) will be added to the oligosaccharide_residue_linkages_.
                std::vector<std::string> oligo_linkages_tokens = gmml::Split(oligosaccharide_linkages_, "\n");
                std::vector<std::string> oligo_name_tokens = gmml::Split(oligosaccharide_name_, "-"); ///According to position of brackets in the name sequence,
                                                                                                      ///the brackets will be added to the oligosaccharide_residue_linkages_
                std::stringstream residue_links_stream;
                std::string full_glycosidic_linkage = "";
                std::vector<std::string> link_tokens = std::vector<std::string>();
                std::string link_left_side = "";
                std::string link_right_side = "";
                size_t end_index = 0;
                std::string first_residue_of_linkage  = "";
                std::string second_residue_of_linkage  = "";

                for(unsigned int i = 0; i < oligo_linkages_tokens.size(); i++) ///Processing linkages line by line
                {
                    full_glycosidic_linkage = oligo_linkages_tokens.at(i);
                    link_tokens = gmml::Split(full_glycosidic_linkage, "}-,");
                    link_left_side = link_tokens.at(1);     ///Getting the first residue of linkage in each line. e.g. {2}RAM(401_A)C1-RAE(402_A)C3, Glycosidic linkage: RAE(402_A)O3
                                                                    ///-> first_mono == RAM(401_A)C1
                    end_index = link_left_side.find_last_of(")");
                    first_residue_of_linkage = link_left_side.substr(0, end_index + 1); ///filtering out atom name. e.g. RAM(401_A)

                    if(oligo_name_tokens.at(i).find("[") != std::string::npos)
                        residue_links_stream << "[" << first_residue_of_linkage;
                    if(oligo_name_tokens.at(i).find("]") != std::string::npos)
                        residue_links_stream << "]-" << first_residue_of_linkage;
                    if(i < oligo_linkages_tokens.size() - 1)
                        residue_links_stream << first_residue_of_linkage << "-";
                    else
                    {
                        link_right_side = link_tokens.at(2);     ///Getting the second residue of linkage in the line. e.g. {1}NAG(1521_A)C1-NAG(1520_A)C4, Glycosidic linkage: NAG(1520_A)O4
                                                                        ///-> second_mono == NAG(1520_A)C4
                        end_index = link_right_side.find_last_of(")");
                        second_residue_of_linkage = link_right_side.substr(0, end_index + 1); ///filtering out atom name. e.g. NAG(1520_A)

                        residue_links_stream << first_residue_of_linkage << "-" << second_residue_of_linkage;
                    }
                }
                oligosaccharide_residue_linkages_ = residue_links_stream.str();
                gmml::FindReplaceString(oligosaccharide_residue_linkages_, "-]", "]");

                out << oligosaccharide_name_;
                gmml::log(__LINE__, __FILE__,  gmml::INF, oligosaccharide_name_);
                out << std::endl;
                out << oligosaccharide_linkages_;
                gmml::log(__LINE__, __FILE__,  gmml::INF, oligosaccharide_linkages_);
                out << std::endl;
//                out << oligosaccharide_residue_linkages_ << std::endl;
            }

            /*! \fn
              * A recursive function in order to generate the name sequence and the linkages between sacchrides of an oligosacchride
              * This functions starts from the root sacchride and generates the name sequence and linkages in reverse order
              * @param oligosaccharide_name The name sequence of the oligosacchride that has been generated by the function so far
              * @param oligosaccharide_linkages The linkages between sacchrides of the oligosacchride which has been generated by the function so far
              * @param i The counter of linkages between the sacchrides of the oligosacchride
              * @param main_root_id The identifier of the main root in the oligosacchride
              * @param is_cycle The boolean variable which represnets whether the sacchrides involved in the oligosacchride forms a cycle
              */
            void GenerateNameLinkage(std::string& oligosaccharide_name, std::string& oligosaccharide_linkages, int& i, int main_root_id, bool& is_cycle)
            {
                std::stringstream name;
                std::stringstream res_linkage;
                std::stringstream res_linkage1;
                if(child_oligos_.size() == 0)
                {
                    if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                        name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
                    else
                        name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
                    oligosaccharide_name = name.str();
                }
                else if(child_oligos_.size() == 1)
                {
                    if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                         name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
                    else
                        name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
                    oligosaccharide_name = name.str();

                    Monosaccharide* mono1 = root_;
                    Oligosaccharide* child_oligo = child_oligos_.at(0);
                    Monosaccharide* mono2 = child_oligo->root_;

                    std::vector<std::string> tokens = gmml::Split(child_oligos_linkages_.at(0), "-");

                    std::stringstream link1;
                    std::string from = tokens.at(0);

                    std::string mono1_ring_atoms_str = mono1->cycle_atoms_str_;
                    std::vector<MolecularModeling::Atom*> mono1_ring_atoms = mono1->cycle_atoms_;
                    std::vector<std::vector<MolecularModeling::Atom*> > mono1_side_atoms = mono1->side_atoms_;
                    int mono1_start_index = 1;
                    std::vector<MolecularModeling::Atom*> mono1_anomeric_side_atoms = mono1_side_atoms.at(0);
                    if(mono1_anomeric_side_atoms.at(0) != NULL)
                        mono1_start_index = 2;
                    if(mono1_ring_atoms_str.find(from) != std::string::npos)
                    {
                        for(std::vector<MolecularModeling::Atom*>::iterator it = mono1_ring_atoms.begin(); it != mono1_ring_atoms.end(); it++)
                        {
                            MolecularModeling::Atom* atom = *it;
                            if(atom != NULL && atom->GetId().compare(from) == 0)
                            {
                                int from_index = std::distance(mono1_ring_atoms.begin(), it) + mono1_start_index;
                                link1 << "-" << gmml::ConvertT(from_index);
                                std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
                                std::vector<std::string> glycosidic_linkage_id_tokens = gmml::Split(tokens.at(1), "_");
                                if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                    res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")"  << atom_id_tokens.at(0);
                                else
                                    res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) <<  "_" << atom_id_tokens.at(3) << ")"  << atom_id_tokens.at(0);

                                if(glycosidic_linkage_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                    res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << ")" << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                else
                                    res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << "_" << glycosidic_linkage_id_tokens.at(3)
                                                 << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                break;
                            }
                        }
                    }
                    else
                    {
                        std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono1_side_atoms.at(mono1_side_atoms.size() - 1);
                        for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
                        {
                            MolecularModeling::Atom* atom = *it;
                            if(atom != NULL && atom->GetId().compare(from) == 0)
                            {
                                int from_index = std::distance(last_carbon_side_atoms.begin(), it) + mono1_side_atoms.size() + mono1_start_index;
                                link1 << "-" << gmml::ConvertT(from_index);
                                std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
                                std::vector<std::string> glycosidic_linkage_id_tokens = gmml::Split(tokens.at(1), "_");
                                if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                    res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")" << atom_id_tokens.at(0);
                                else
                                    res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);

                                if(glycosidic_linkage_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                    res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << ")" << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                else
                                    res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << "_" << glycosidic_linkage_id_tokens.at(3)
                                                 << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                break;
                            }
                        }
                    }

                    std::stringstream link;
                    std::string to = tokens.at(2);

                    std::string mono2_ring_atoms_str = mono2->cycle_atoms_str_;
                    std::vector<MolecularModeling::Atom*> mono2_ring_atoms = mono2->cycle_atoms_;
                    std::vector<std::vector<MolecularModeling::Atom*> > mono2_side_atoms = mono2->side_atoms_;
                    std::vector<MolecularModeling::Atom*> mono2_anomeric_side_atoms = mono2_side_atoms.at(0);

                    int mono2_start_index = 1;
                    if(mono2_anomeric_side_atoms.at(0) != NULL)
                        mono2_start_index = 2;
                    if(mono2_ring_atoms_str.find(to) != std::string::npos)
                    {
                        for(std::vector<MolecularModeling::Atom*>::iterator it = mono2_ring_atoms.begin(); it != mono2_ring_atoms.end(); it++)
                        {
                            MolecularModeling::Atom* atom = *it;
                            if(atom != NULL && atom->GetId().compare(to) == 0)
                            {
                                int to_index = std::distance(mono2_ring_atoms.begin(), it) + mono2_start_index;
                                link << gmml::ConvertT(to_index);
                                std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
                                if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                    res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")" << atom_id_tokens.at(0);
                                else
                                    res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);
                                break;
                            }
                        }
                    }
                    else
                    {
                        std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono2_side_atoms.at(mono2_side_atoms.size() - 1);
                        for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
                        {
                            MolecularModeling::Atom* atom = *it;
                            if(atom != NULL && atom->GetId().compare(to) == 0)
                            {
                                int to_index = std::distance(last_carbon_side_atoms.begin(), it) + mono2_side_atoms.size() + mono2_start_index;
                                link << gmml::ConvertT(to_index);
                                std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
                                if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                    res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")" << atom_id_tokens.at(0);
                                else
                                    res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);
                                break;
                            }
                        }
                    }

                    link << link1.str() << oligosaccharide_name;
                    oligosaccharide_name = link.str();
                    res_linkage << res_linkage1.str() << oligosaccharide_linkages;
                    oligosaccharide_linkages = res_linkage.str();

                    i++;
                    if(main_root_id == child_oligo->root_->mono_id)
                        is_cycle = true;
                    child_oligo->GenerateNameLinkage(oligosaccharide_name, oligosaccharide_linkages, i, main_root_id, is_cycle);
                }
                else
                {
                    std::map<int, int> root_child_index_map = std::map<int, int>();

                    int mono1_start_index = 1;
                    if(root_->side_atoms_.at(0).at(0) != NULL)///mono has a carbon at position -1 so the ring atom indices should start from 2
                        mono1_start_index = 2;
                    for(std::vector<std::string>::iterator it = child_oligos_linkages_.begin(); it != child_oligos_linkages_.end(); it++)
                    {
                        std::string child_linkage = (*it);
                        int linkage_index = std::distance(child_oligos_linkages_.begin(), it);
                        std::vector<std::string> child_linkage_tokens = gmml::Split(child_linkage, "-");
                        std::string from = child_linkage_tokens.at(0);///the atom id from current mono that is involved in the linkage

                        if(root_->cycle_atoms_str_.find(from) != std::string::npos)///if the atom that is involved in the linkage is one of the current mono cycle atoms
                        {
                            for(std::vector<MolecularModeling::Atom*>::iterator it = root_->cycle_atoms_.begin(); it != root_->cycle_atoms_.end(); it++)
                            {
                                MolecularModeling::Atom* atom = *it;
                                if(atom != NULL && atom->GetId().compare(from) == 0)
                                {
                                    int from_index = std::distance(root_->cycle_atoms_.begin(), it) + mono1_start_index;///index position of the atom in the mono
                                    root_child_index_map[linkage_index] = from_index;
                                    break;
                                }
                            }
                        }
                        else///the atom that is involved in the linkage is one of the current mono side atoms
                        {
                            std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = root_->side_atoms_.at(root_->side_atoms_.size() - 1);
                            for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
                            {
                                MolecularModeling::Atom* atom = *it;
                                if(atom != NULL && atom->GetId().compare(from) == 0)
                                {
                                    int from_index = std::distance(last_carbon_side_atoms.begin(), it) + root_->side_atoms_.size() + mono1_start_index;///index position of the side atom in the mono
                                    root_child_index_map[linkage_index] = from_index;
                                    break;
                                }
                            }
                        }
                    }

                    bool isFirstChild = true;
                    while(root_child_index_map.size() > 0)
                    {
                        res_linkage.str("");
                        res_linkage1.str("");
                        int min = 100;
                        int index_min = 0;///index of the child oligo with the min children in child_oligos
                        for(std::map<int, int>::iterator it0 = root_child_index_map.begin(); it0 != root_child_index_map.end(); it0++)
                        {
                            int map_key = (*it0).first;
                            int map_value = (*it0).second;
                            if(map_value < min)
                            {
                                index_min = map_key;
                                min = map_value;
                            }
                        }
                        if(isFirstChild)
                        {
                            isFirstChild = false;
                            if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                                name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
                            else
                                name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
                            oligosaccharide_name = name.str();
                        }

                        Monosaccharide* mono1 = root_;
                        Oligosaccharide* child_oligo = child_oligos_.at(index_min);
                        Monosaccharide* mono2 = child_oligo->root_;

                        std::vector<std::string> tokens = gmml::Split(child_oligos_linkages_.at(index_min), "-");

                        std::stringstream link1;
                        std::string from = tokens.at(0);

                        std::string mono1_ring_atoms_str = mono1->cycle_atoms_str_;
                        std::vector<MolecularModeling::Atom*> mono1_ring_atoms = mono1->cycle_atoms_;
                        std::vector<std::vector<MolecularModeling::Atom*> > mono1_side_atoms = mono1->side_atoms_;
                        int mono1_start_index = 1;
                        std::vector<MolecularModeling::Atom*> mono1_anomeric_side_atoms = mono1_side_atoms.at(0);
                        if(mono1_anomeric_side_atoms.at(0) != NULL)
                            mono1_start_index = 2;
                        if(mono1_ring_atoms_str.find(from) != std::string::npos)
                        {
                            for(std::vector<MolecularModeling::Atom*>::iterator it = mono1_ring_atoms.begin(); it != mono1_ring_atoms.end(); it++)
                            {
                                MolecularModeling::Atom* atom = *it;
                                if(atom != NULL && atom->GetId().compare(from) == 0)
                                {
                                    int from_index = std::distance(mono1_ring_atoms.begin(), it) + mono1_start_index;
                                    link1 << "-" << gmml::ConvertT(from_index);
                                    std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
                                    std::vector<std::string> glycosidic_linkage_id_tokens = gmml::Split(tokens.at(1), "_");
                                    if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                        res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")"  << atom_id_tokens.at(0);
                                    else
                                        res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);

                                    if(glycosidic_linkage_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                        res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << ")" << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                    else
                                        res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << "_" << glycosidic_linkage_id_tokens.at(3)
                                                     << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono1_side_atoms.at(mono1_side_atoms.size() - 1);
                            for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
                            {
                                MolecularModeling::Atom* atom = *it;
                                if(atom != NULL && atom->GetId().compare(from) == 0)
                                {
                                    int from_index = std::distance(last_carbon_side_atoms.begin(), it) + mono1_side_atoms.size() + mono1_start_index;
                                    link1 << "-" << gmml::ConvertT(from_index);
                                    std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
                                    std::vector<std::string> glycosidic_linkage_id_tokens = gmml::Split(tokens.at(1), "_");
                                    if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                        res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")" << atom_id_tokens.at(0);
                                    else
                                        res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);

                                    if(glycosidic_linkage_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                        res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << ")" << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                    else
                                        res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << "_" << glycosidic_linkage_id_tokens.at(3)
                                                     << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                    break;
                                }
                            }
                        }

                        std::stringstream link;
                        std::string to = tokens.at(2);

                        std::string mono2_ring_atoms_str = mono2->cycle_atoms_str_;
                        std::vector<MolecularModeling::Atom*> mono2_ring_atoms = mono2->cycle_atoms_;
                        std::vector<std::vector<MolecularModeling::Atom*> > mono2_side_atoms = mono2->side_atoms_;
                        std::vector<MolecularModeling::Atom*> mono2_anomeric_side_atoms = mono2_side_atoms.at(0);

                        int mono2_start_index = 1;
                        if(mono2_anomeric_side_atoms.at(0) != NULL)
                            mono2_start_index = 2;
                        if(mono2_ring_atoms_str.find(to) != std::string::npos)
                        {
                            for(std::vector<MolecularModeling::Atom*>::iterator it = mono2_ring_atoms.begin(); it != mono2_ring_atoms.end(); it++)
                            {
                                MolecularModeling::Atom* atom = *it;
                                if(atom != NULL && atom->GetId().compare(to) == 0)
                                {
                                    int to_index = std::distance(mono2_ring_atoms.begin(), it) + mono2_start_index;
                                    link << gmml::ConvertT(to_index);
                                    std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
                                    if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                        res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")" << atom_id_tokens.at(0);
                                    else
                                        res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);
                                    break;
                                }
                            }
                        }
                        else
                        {
                            std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono2_side_atoms.at(mono2_side_atoms.size() - 1);
                            for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
                            {
                                MolecularModeling::Atom* atom = *it;
                                if(atom != NULL && atom->GetId().compare(to) == 0)
                                {
                                    int to_index = std::distance(last_carbon_side_atoms.begin(), it) + mono2_side_atoms.size() + mono2_start_index;
                                    link << gmml::ConvertT(to_index);
                                    std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
                                    if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
                                        res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")" << atom_id_tokens.at(0);
                                    else
                                        res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);
                                    break;
                                }
                            }
                        }

                        if(root_child_index_map.size() > 1)
                            link << link1.str() << "]" << oligosaccharide_name;
                        else
                            link << link1.str() << oligosaccharide_name;
                        oligosaccharide_name = link.str();
                        res_linkage << res_linkage1.str() << oligosaccharide_linkages;
                        oligosaccharide_linkages = res_linkage.str();

                        i++;
                        child_oligo->GenerateNameLinkage(oligosaccharide_name, oligosaccharide_linkages, i, main_root_id, is_cycle);
                        std::stringstream name;
                        if(root_child_index_map.size() > 1)
                            name << "[" << oligosaccharide_name;
                        else
                            name << oligosaccharide_name;
                        oligosaccharide_name = name.str();
                        root_child_index_map.erase(index_min);
                    }
                }
            }
    };
}

#endif // OLIGOSACCHARIDE_HPP
