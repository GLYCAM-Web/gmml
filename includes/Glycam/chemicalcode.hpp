#ifndef CHEMICALCODE_HPP
#define CHEMICALCODE_HPP

#include <string>
#include <vector>
#include <sstream>

#include "../MolecularModeling/atom.hpp"

namespace Glycam
{
    struct ChemicalCode {
            std::string base_;
            std::vector<std::string> left_up_;
            std::vector<MolecularModeling::Atom*> left_up_atoms_;
            std::vector<std::string> left_down_;
            std::vector<MolecularModeling::Atom*> left_down_atoms_;
            std::vector<std::string> left_middle_;
            std::vector<MolecularModeling::Atom*> left_middle_atoms_;
            std::vector<std::string> right_middle_;
            std::vector<MolecularModeling::Atom*> right_middle_atoms_;
            std::vector<std::string> right_up_;
            std::vector<MolecularModeling::Atom*> right_up_atoms_;
            std::vector<std::string> right_down_;
            std::vector<MolecularModeling::Atom*> right_down_atoms_;

            std::string toString()
            {
                std::stringstream ss;
                for(std::vector<std::string>::iterator it = this->left_down_.begin(); it != this->left_down_.end(); it++)
                {
                    ss << "_" << (*it);
                }
                for(std::vector<std::string>::iterator it = this->left_up_.begin(); it != this->left_up_.end(); it++)
                {
                    ss << "^" << (*it);
                }
                for(std::vector<std::string>::iterator it = this->left_middle_.begin(); it != this->left_middle_.end(); it++)
                {
                    ss << (*it);
                }
                ss << this->base_;
                for(std::vector<std::string>::iterator it = this->right_middle_.begin(); it != this->right_middle_.end(); it++)
                {
                    ss << (*it);
                }
                for(std::vector<std::string>::iterator it = this->right_down_.begin(); it != this->right_down_.end(); it++)
                {
                    ss << "_" << (*it);
                }
                for(std::vector<std::string>::iterator it = this->right_up_.begin(); it != this->right_up_.end(); it++)
                {
                    ss << "^" << (*it);
                }
                return ss.str();
            }

            void Print(std::ostream& out)
            {
                out << std::setw(8) << std::right;
                std::stringstream lu;
                for(int i = 0; i < this->left_up_.size(); i++)
                {
                    lu << this->left_up_.at(i) << " ";
                }
                out << lu.str();

                out << std::setw(8) << std::right;
                std::stringstream ru;
                for(int i = 0; i < this->right_up_.size(); i++)
                {
                    ru << this->right_up_.at(i) << " ";
                }
                out << ru.str();
                out << std::endl;

                out << std::setw(6) << std::right;
                std::stringstream lm;
                for(int i = 0; i < this->left_middle_.size(); i++)
                {
                    lm << this->left_middle_.at(i) << " ";
                }
                out << lm.str();

                out << std::setw(3) << std::right <<  this->base_;

                out << std::setw(6) << std::right;
                std::stringstream rm;
                for(int i = 0; i < this->right_middle_.size(); i++)
                {
                    rm << this->right_middle_.at(i) << " ";
                }
                out << rm.str();
                out << std::endl;

                out << std::setw(8) << std::right;
                std::stringstream ld;
                for(int i = 0; i < this->left_down_.size(); i++)
                {
                    ld << this->left_down_.at(i) << " ";
                }
                out << ld.str();

                out << std::setw(8) << std::right;
                std::stringstream rd;
                for(int i = 0; i < this->right_down_.size(); i++)
                {
                    rd << this->right_down_.at(i) << " ";
                }
                out << rd.str();
                out << std::endl;
            }
    } ;
}

#endif // CHEMICALCODE_HPP
