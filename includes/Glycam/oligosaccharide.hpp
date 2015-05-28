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
            std::vector<std::string> linkages_;
            std::vector<Monosaccharide*> monos_;
            std::string oligosaccharide_name_;
            std::string residue_linkages_;

            void Print(std::ostream& out = std::cout)
            {
                std::stringstream name;
                std::stringstream res_linkage;
                for(int i = 0; i < monos_.size() - 1; i++)
                {
                    Monosaccharide* mono1 = monos_.at(i);
                    Monosaccharide* mono2 = monos_.at(i+1);
                    if(mono1->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                        name << mono1->sugar_name_.monosaccharide_short_name_;
                    else if(mono1->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                        name << mono1->sugar_name_.monosaccharide_stereochemistry_short_name_;


                    std::vector<std::string> tokens = gmml::Split(linkages_.at(i), "-");
                    std::string from = tokens.at(0);

                    std::string mono1_ring_atoms_str = mono1->cycle_atoms_str_;
                    std::vector<MolecularModeling::Atom*> mono1_ring_atoms = mono1->cycle_atoms_;
                    std::vector<std::vector<MolecularModeling::Atom*> > mono1_side_atoms = mono1->side_atoms_;

                    if(mono1_ring_atoms_str.find(from) != std::string::npos)
                    {
                        for(std::vector<MolecularModeling::Atom*>::iterator it = mono1_ring_atoms.begin(); it != mono1_ring_atoms.end(); it++)
                        {
                            MolecularModeling::Atom* atom = *it;
                            if(atom != NULL && atom->GetId().compare(from) == 0)
                            {
                                int from_index = std::distance(mono1_ring_atoms.begin(), it) + 1;
                                name << gmml::ConvertT(from_index) << "-";
                                std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
                                res_linkage << "[" << i+1 << "]" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")"  << atom_id_tokens.at(0) << "-";
                                break;
                            }
                        }
                    }
                    else
                    {
                        std::vector<MolecularModeling::Atom*> anomeric_side_atoms = mono1_side_atoms.at(0);
                        if(anomeric_side_atoms.at(0) != NULL && anomeric_side_atoms.at(0)->GetId().compare(from) == 0)
                            name << "0-";
                        else
                        {
                            std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono1_side_atoms.at(mono1_side_atoms.size() - 1);
                            for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
                            {
                                MolecularModeling::Atom* atom = *it;
                                if(atom != NULL && atom->GetId().compare(from) == 0)
                                {
                                    int from_index = std::distance(last_carbon_side_atoms.begin(), it) + mono1_side_atoms.size() + 1;
                                    name << gmml::ConvertT(from_index) << "-";
                                    std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
                                    res_linkage << "[" << i+1 << "]" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0) << "-";
                                    break;
                                }
                            }

                        }
                    }

                    std::string to = tokens.at(2);
                    std::string mono2_ring_atoms_str = mono2->cycle_atoms_str_;
                    std::vector<MolecularModeling::Atom*> mono2_ring_atoms = mono2->cycle_atoms_;
                    std::vector<std::vector<MolecularModeling::Atom*> > mono2_side_atoms = mono2->side_atoms_;
                    if(mono2_ring_atoms_str.find(to) != std::string::npos)
                    {
                        for(std::vector<MolecularModeling::Atom*>::iterator it = mono2_ring_atoms.begin(); it != mono2_ring_atoms.end(); it++)
                        {
                            MolecularModeling::Atom* atom = *it;
                            if(atom != NULL && atom->GetId().compare(to) == 0)
                            {
                                int to_index = std::distance(mono2_ring_atoms.begin(), it) + 1;
                                name << gmml::ConvertT(to_index);
                                std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
                                res_linkage << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0) << std::endl;
                                break;
                            }
                        }
                    }
                    else
                    {
                        std::vector<MolecularModeling::Atom*> anomeric_side_atoms = mono2_side_atoms.at(0);
                        if(anomeric_side_atoms.at(0) != NULL && anomeric_side_atoms.at(0)->GetId().compare(to) == 0)
                            name << "0";
                        else
                        {
                            std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono2_side_atoms.at(mono2_side_atoms.size() - 1);
                            for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
                            {
                                MolecularModeling::Atom* atom = *it;
                                if(atom != NULL && atom->GetId().compare(to) == 0)
                                {
                                    int to_index = std::distance(last_carbon_side_atoms.begin(), it) + mono2_side_atoms.size() + 1;
                                    name << gmml::ConvertT(to_index);
                                    std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
                                    res_linkage << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << atom_id_tokens.at(3) << ")" << atom_id_tokens.at(0) << std::endl;
                                    break;
                                }
                            }
                        }
                    }

//                    name << "(" << linkages_.at(i) << ")";

                }
                Monosaccharide* mono1 = monos_.at(monos_.size() - 1);
                Monosaccharide* mono2 = monos_.at(0);
                if(mono1->cycle_atoms_str_.compare(mono2->cycle_atoms_str_) != 0)
                {
                    if(mono1->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                        name << mono1->sugar_name_.monosaccharide_short_name_;
                    else if(mono1->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
                        name << mono1->sugar_name_.monosaccharide_stereochemistry_short_name_;
                }
                oligosaccharide_name_ = name.str();
                residue_linkages_ = res_linkage.str();

                out << name.str() << std::endl;
                out << res_linkage.str() << std::endl;
                out << std::endl;
            }
    };
}

#endif // OLIGOSACCHARIDE_HPP
