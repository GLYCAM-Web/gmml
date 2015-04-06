#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

#include <string>
#include <vector>
#include <sstream>

namespace gmml
{
    struct SugarName {
      std::string name_;
      std::string isomer_;
      std::string configuration_;
      std::string ring_type_;
      std::string monosaccharide_name_;
    } ;
    struct ChemicalCode {
            std::string base_;
            std::vector<std::string> left_up_;
            std::vector<std::string> left_down_;
            std::vector<std::string> right_up_;
            std::vector<std::string> right_down_;

            void Print(std::ostream& out)
            {

                out << std::setw(6) << std::right;
                std::stringstream lu;
                for(int i = 0; i < this->left_up_.size(); i++)
                {
                    lu << this->left_up_.at(i) << " ";
                }
                out << lu.str();

                out << std::setw(6) << std::right;
                std::stringstream ru;
                for(int i = 0; i < this->right_up_.size(); i++)
                {
                    ru << this->right_up_.at(i) << " ";
                }
                out << ru.str();
                out << std::endl;

                out << std::setw(7) << std::right <<  this->base_ << std::endl;

                out << std::setw(6) << std::right;
                std::stringstream ld;
                for(int i = 0; i < this->left_down_.size(); i++)
                {
                    ld << this->left_down_.at(i) << " ";
                }
                out << ld.str();

                out << std::setw(6) << std::right;
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

#endif // STRUCTURES_HPP
