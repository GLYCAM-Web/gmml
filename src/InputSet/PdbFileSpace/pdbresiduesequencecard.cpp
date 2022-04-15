// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbResidueSequenceCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueSequenceCard::PdbResidueSequenceCard() : chain_id_(' '), number_of_residues_(gmml::dNotSet) {}

PdbResidueSequenceCard::PdbResidueSequenceCard(char chain_id, int number_of_residues, const std::vector<std::string> &residue_names) : chain_id_(chain_id),
    number_of_residues_(number_of_residues)
{
    residue_names_.clear();
    for(std::vector<std::string>::const_iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
    }
}

PdbResidueSequenceCard::PdbResidueSequenceCard(std::stringstream& stream_block)
{
    std::string line;
    bool is_chain_id_set = false, is_number_of_residues_set = false;
    std::stringstream ss;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_chain_id_set){
            if(line.substr(11,1) == " ")
                chain_id_ = ' ';
            else
                chain_id_ = gmml::ConvertString<char>(line.substr(11,1));
            is_chain_id_set=true;
        }
        if(!is_number_of_residues_set){
            if(line.substr(13,4) == "    ")
                number_of_residues_ = gmml::iNotSet;
            else
                number_of_residues_ = gmml::ConvertString<int>(line.substr(13,4));
            is_number_of_residues_set=true;
        }
        ss << line.substr(19,51) << " ";

        getline(stream_block, line);
        temp = line;
    }
    residue_names_ = gmml::Split(ss.str(), " ");
    for(std::vector<std::string>::iterator it = residue_names_.begin(); it != residue_names_.end(); it++)
    {
        gmml::Trim(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbResidueSequenceCard::GetChainId()
{
    return chain_id_;
}

int PdbResidueSequenceCard::GetNumberOfResidues()
{
    return number_of_residues_;
}

std::vector<std::string> PdbResidueSequenceCard::GetResidueNames()
{
    return residue_names_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueSequenceCard::SetChainId(char chain_id)
{
    chain_id_ = chain_id;
}

void PdbResidueSequenceCard::SetNumberOfResidues(int number_of_residues)
{
    number_of_residues_ = number_of_residues;
}

void PdbResidueSequenceCard::SetResidueNames(const std::vector<std::string> residue_names)
{
    residue_names_.clear();
    for(std::vector<std::string>::const_iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
    }
}

void PdbResidueSequenceCard::AddResidueName(const std::string residue_name)
{
    residue_names_.push_back(residue_name);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueSequenceCard::Print(std::ostream &out)
{
    out << "Chain ID: " << chain_id_
        << ", " << "Number of Residues: ";
    if(number_of_residues_ != gmml::iNotSet)
        out << number_of_residues_;
    else
        out << " ";
    out << "Residue Names: ";
    for(std::vector<std::string>::iterator it = residue_names_.begin(); it != residue_names_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << std::endl;
}
