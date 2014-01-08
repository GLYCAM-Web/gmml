// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbresiduesequence.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueSequence::PdbResidueSequence() : chain_id_(' '), number_of_residues_(kNotSet) {}

PdbResidueSequence::PdbResidueSequence(char chain_id, int number_of_residues, const vector<string> &residue_names) : chain_id_(chain_id),
    number_of_residues_(number_of_residues)
{
    residue_names_.clear();
    for(vector<string>::const_iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
    }
}

PdbResidueSequence::PdbResidueSequence(stringstream& stream_block)
{
    string line;
    bool is_chain_id_set = false, is_number_of_residues_set = false;
    stringstream ss;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_chain_id_set){
            chain_id_ = ConvertString<char>(line.substr(11,1));
            is_chain_id_set=true;
        }
        if(!is_number_of_residues_set){
            number_of_residues_ = ConvertString<int>(line.substr(13,4));
            is_number_of_residues_set=true;
        }
        ss << line.substr(19,51) << " ";

        getline(stream_block, line);
    }
    residue_names_ = Split(ss.str(), " ");
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbResidueSequence::GetChainId()
{
    return chain_id_;
}

int PdbResidueSequence::GetNumberOfResidues()
{
    return number_of_residues_;
}

vector<string> PdbResidueSequence::GetResidueNames()
{
    return residue_names_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueSequence::SetChainId(char chain_id)
{
    chain_id_ = chain_id;
}

void PdbResidueSequence::SetNumberOfResidues(int number_of_residues)
{
    number_of_residues_ = number_of_residues;
}

void PdbResidueSequence::SetResidueNames(const vector<string> residue_names)
{
    residue_names_.clear();
    for(vector<string>::const_iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
    }
}

void PdbResidueSequence::AddResidueName(const string residue_name)
{
    residue_names_.push_back(residue_name);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

