#ifndef OLIGOSACCHARIDE_HPP
#define OLIGOSACCHARIDE_HPP

#include <vector>
#include <string>
#include <sstream>

#include "../MolecularModeling/assembly.hpp"
#include "../utils.hpp"

namespace Glycam
{
    struct Oligosaccharide
    {
            Monosaccharide* root_;
            std::vector<Oligosaccharide*> child_oligos_;
            std::vector<std::string> child_oligos_linkages_;
            std::string oligosaccharide_name_;
            std::string oligosaccharide_linkages_;

            void Print(std::string terminal = "", std::ostream& out = std::cout)
            {
                oligosaccharide_name_ = "";
                oligosaccharide_linkages_ = "";
                int linkage_index = 0;
                std::stringstream ss;
                if(terminal.compare("") != 0)
                    ss << "-" << terminal;
                oligosaccharide_name_ = ss.str();
//                std::cout << root_->mono_id << std::endl;
                bool is_cycle = false;
                this->GenerateNameLinkage(oligosaccharide_name_, oligosaccharide_linkages_, linkage_index, root_->mono_id, is_cycle);

                if(is_cycle)
                {
                    int start_index = oligosaccharide_name_.find_first_of('-') + 1;
                    int end_index = oligosaccharide_name_.find_last_of('-');
                    std::vector<std::string> tokens = gmml::Split(oligosaccharide_name_, "-");

                    std::string middle_sub_name = oligosaccharide_name_.substr(start_index, end_index - start_index);
                    std::stringstream new_name;
                    new_name << tokens.at(0) << "-[" << middle_sub_name << "]-" << tokens.at(tokens.size() - 1).at(0);
                    oligosaccharide_name_ = new_name.str();
                }
                out << oligosaccharide_name_ << std::endl;
                out << oligosaccharide_linkages_ << std::endl;
            }

            void GenerateNameLinkage(std::string& oligosaccharide_name, std::string& oligosaccharide_linkages, int& i, int main_root_id, bool& is_cycle)
            {
                std::stringstream name;
                std::stringstream res_linkage;
                std::stringstream res_linkage1;
                if(child_oligos_.size() == 0)
                {
//                    std::cout << root_->mono_id << " has 0 child" << std::endl;
                    if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                        name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
                    else
                        name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
                    oligosaccharide_name = name.str();
//                    std::cout << "NAME: " << oligosaccharide_name << std::endl;
                }
                else if(child_oligos_.size() == 1)
                {
//                    std::cout << root_->mono_id << " has 1 child" << std::endl;
                    if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                         name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
                    else
                        name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
                    oligosaccharide_name = name.str();
//                    std::cout << "NAME: " << oligosaccharide_name << std::endl;

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
                                    res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")"  << atom_id_tokens.at(0)
                                                 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << glycosidic_linkage_id_tokens.at(3)
                                                 << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                else
                                    res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")"  << atom_id_tokens.at(0)
                                                 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << glycosidic_linkage_id_tokens.at(3)
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
                                    res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")" << atom_id_tokens.at(0)
                                                 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << glycosidic_linkage_id_tokens.at(3)
                                                 << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                else
                                    res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0)
                                                 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << glycosidic_linkage_id_tokens.at(3)
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
                                    res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);
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
                                    res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);
                                break;
                            }
                        }
                    }

                    link << link1.str() << oligosaccharide_name;
                    oligosaccharide_name = link.str();
//                    std::cout << "NAME: " << oligosaccharide_name << std::endl;
                    res_linkage << res_linkage1.str() << oligosaccharide_linkages;
                    oligosaccharide_linkages = res_linkage.str();

                    i++;
//                    std::cout << child_oligo->root_->mono_id << std::endl;
                    if(main_root_id == child_oligo->root_->mono_id)
                        is_cycle = true;
                    child_oligo->GenerateNameLinkage(oligosaccharide_name, oligosaccharide_linkages, i, main_root_id, is_cycle);
                }
                else
                {
//                    std::cout << root_->mono_id << " has more than 1 child" << std::endl;
                    std::map<int, int> root_child_index_map = std::map<int, int>();

                    int mono1_start_index = 1;
                    if(root_->side_atoms_.at(0).at(0) != NULL)///mono has a carbon at position -1 so the ring atom indices should start from 2
                        mono1_start_index = 2;
                    for(std::vector<std::string>::iterator it = child_oligos_linkages_.begin(); it != child_oligos_linkages_.end(); it++)
                    {
                        std::string child_linkage = (*it);
//                        int linkage_index = std::distance(child_oligos_linkages_.begin(), it) + child_oligos_linkages_.size() ;
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
//                        std::cout << "mono" << root_->mono_id << " root_child_map size " << root_child_index_map.size() << std::endl;
                        int min = 100;
                        int index_min = 0;///index of the child oligo with the min children in child_oligos
//                        std::cout << "root_child_index_map" << std::endl;
                        for(std::map<int, int>::iterator it0 = root_child_index_map.begin(); it0 != root_child_index_map.end(); it0++)
                        {
                            int map_key = (*it0).first;
                            int map_value = (*it0).second;
//                            std::cout << "key: " << map_key << ", value: " << map_value << std::endl;
                            if(map_value < min)
                            {
                                index_min = map_key;
                                min = map_value;
                            }
                        }

//                        std::cout << "chosen child mono is " << child_oligos_.at(index_min)->root_->mono_id << std::endl;

                        if(isFirstChild)
                        {
                            isFirstChild = false;
                            if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                                name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
                            else
                                name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
                            oligosaccharide_name = name.str();
                        }
//                        std::cout << "NAME: " << oligosaccharide_name << std::endl;

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
                                        res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")"  << atom_id_tokens.at(0)
                                                     <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << glycosidic_linkage_id_tokens.at(3)
                                                     << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                    else
                                        res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0)
                                                     <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << glycosidic_linkage_id_tokens.at(3)
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
                                        res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << ")" << atom_id_tokens.at(0)
                                                     <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << glycosidic_linkage_id_tokens.at(3)
                                                     << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
                                    else
                                        res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0)
                                                     <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << glycosidic_linkage_id_tokens.at(3)
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
                                        res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);
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
                                        res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0);
                                    break;
                                }
                            }
                        }

                        if(root_child_index_map.size() > 1)
                            link << link1.str() << "]" << oligosaccharide_name;
                        else
                            link << link1.str() << oligosaccharide_name;
                        oligosaccharide_name = link.str();
//                        std::cout << "NAME: " << oligosaccharide_name << std::endl;
                        res_linkage << res_linkage1.str() << oligosaccharide_linkages;
                        oligosaccharide_linkages = res_linkage.str();

                        i++;
//                        std::cout << "call function " <<std::endl;
//                        std::cout << child_oligo->root_->mono_id << std::endl;
                        child_oligo->GenerateNameLinkage(oligosaccharide_name, oligosaccharide_linkages, i, main_root_id, is_cycle);
                        std::stringstream name;
                        if(root_child_index_map.size() > 1)
                            name << "[" << oligosaccharide_name;
                        else
                            name << oligosaccharide_name;
                        oligosaccharide_name = name.str();
//                        std::cout << "NAME: " << oligosaccharide_name << std::endl;
//                        std::cout << "earasing index min " << index_min << " mono" << child_oligos_.at(index_min)->root_->mono_id << std::endl;
                        root_child_index_map.erase(index_min);

                    }
                }
            }


    };
}

#endif // OLIGOSACCHARIDE_HPP
