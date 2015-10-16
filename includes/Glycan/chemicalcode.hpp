#ifndef CHEMICALCODE_HPP
#define CHEMICALCODE_HPP

#include <string>
#include <vector>
#include <sstream>

#include "../MolecularModeling/atom.hpp"

namespace Glycan
{
    struct ChemicalCode {
            std::string base_;                                              /*!< Base of the chemical code which can be P/F which represnets pyranose/furanose accordingly >*/
            std::vector<std::string> left_up_;                              /*!< The list of ring indices of the atoms that are at the left side of the ring and their position is UP with respect to the ring >*/
            std::vector<MolecularModeling::Atom*> left_up_atoms_;           /*!< The list of atoms that are at the left side of the ring and their position is UP with respect to the ring >*/
            std::vector<std::string> left_down_;                            /*!< The list of ring indices of the atoms that are at the left side of the ring and their position is DOWN with respect to the ring >*/
            std::vector<MolecularModeling::Atom*> left_down_atoms_;         /*!< The list of atoms that are at the left side of the ring and their position is DOWN with respect to the ring >*/
            std::vector<std::string> left_middle_;                          /*!< The list of ring indices of the atoms that are at the left side of the ring and are deoxy >*/
            std::vector<MolecularModeling::Atom*> left_middle_atoms_;       /*!< The list of atoms that are at the left side of the ring and are deoxy >*/
            std::vector<std::string> right_middle_;                         /*!< The list of ring indices of the atoms that are at the right side of the ring and are deoxy >*/
            std::vector<MolecularModeling::Atom*> right_middle_atoms_;      /*!< The list of atoms that are at the right side of the ring and are deoxy >*/
            std::vector<std::string> right_up_;                             /*!< The list of ring indices of the atoms that are at the right side of the ring and their position is UP with respect to the ring >*/
            std::vector<MolecularModeling::Atom*> right_up_atoms_;          /*!< The list of atoms that are at the right side of the ring and their position is UP with respect to the ring >*/
            std::vector<std::string> right_down_;                           /*!< The list of ring indices of the atoms that are at the right side of the ring and their position is DOWN with respect to the ring >*/
            std::vector<MolecularModeling::Atom*> right_down_atoms_;        /*!< The list of atoms that are at the right side of the ring and their position is DOWN with respect to the ring >*/

            /*! \fn
              * A function in order to convert the chemical code structure into the string version of it
              * @return ss The string version of the chemical code structure
              */
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

            /*! \fn
              * A function to print out the chemical code in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out)
            {
                out << std::setw(8) << std::right;
                std::stringstream lu;
                for(unsigned int i = 0; i < this->left_up_.size(); i++)
                {
                    lu << this->left_up_.at(i) << " ";
                }
                out << lu.str();

                out << std::setw(8) << std::right;
                std::stringstream ru;
                for(unsigned int i = 0; i < this->right_up_.size(); i++)
                {
                    ru << this->right_up_.at(i) << " ";
                }
                out << ru.str();
                out << std::endl;

                out << std::setw(6) << std::right;
                std::stringstream lm;
                for(unsigned int i = 0; i < this->left_middle_.size(); i++)
                {
                    lm << this->left_middle_.at(i) << " ";
                }
                out << lm.str();

                out << std::setw(3) << std::right <<  this->base_;

                out << std::setw(6) << std::right;
                std::stringstream rm;
                for(unsigned int i = 0; i < this->right_middle_.size(); i++)
                {
                    rm << this->right_middle_.at(i) << " ";
                }
                out << rm.str();
                out << std::endl;

                out << std::setw(8) << std::right;
                std::stringstream ld;
                for(unsigned int i = 0; i < this->left_down_.size(); i++)
                {
                    ld << this->left_down_.at(i) << " ";
                }
                out << ld.str();

                out << std::setw(8) << std::right;
                std::stringstream rd;
                for(unsigned int i = 0; i < this->right_down_.size(); i++)
                {
                    rd << this->right_down_.at(i) << " ";
                }
                out << rd.str();
                out << std::endl;
            }
    } ;
}

#endif // CHEMICALCODE_HPP
