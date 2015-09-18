// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequence.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueSequence::PdbResidueSequence() : chain_id_(' '), number_of_residues_(dNotSet) {}

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
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_chain_id_set){
            if(line.substr(11,1) == " ")
                chain_id_ = ' ';
            else
                chain_id_ = ConvertString<char>(line.substr(11,1));
            is_chain_id_set=true;
        }
        if(!is_number_of_residues_set){
            if(line.substr(13,4) == "    ")
                number_of_residues_ = iNotSet;
            else
                number_of_residues_ = ConvertString<int>(line.substr(13,4));
            is_number_of_residues_set=true;
        }
        ss << line.substr(19,51) << " ";

        getline(stream_block, line);
        temp = line;
    }
    residue_names_ = Split(ss.str(), " ");
    for(vector<string>::iterator it = residue_names_.begin(); it != residue_names_.end(); it++)
    {
        Trim(*it);
    }
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
void PdbResidueSequence::Print(ostream &out)
{
    out << "Chain ID: " << chain_id_
        << ", " << "Number of Residues: ";
    if(number_of_residues_ != iNotSet)
        out << number_of_residues_;
    else
        out << " ";
    out << "Residue Names: ";
    for(vector<string>::iterator it = residue_names_.begin(); it != residue_names_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << endl;
}
