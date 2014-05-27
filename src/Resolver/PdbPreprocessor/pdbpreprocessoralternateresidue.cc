
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorAlternateResidue::PdbPreprocessorAlternateResidue() {}

PdbPreprocessorAlternateResidue::PdbPreprocessorAlternateResidue(string residue_name, char chain_id, int sequence_number, char residue_insertion_code,
                                                                 vector<char> residue_alternate_location, vector<bool> selected_alternate_location) :
    residue_name_(residue_name), residue_chain_id_(chain_id), residue_sequence_number_(sequence_number), residue_insertion_code_(residue_insertion_code)
{
    residue_alternate_location_ = vector<char>();
    for(vector<char>::iterator it = residue_alternate_location.begin(); it != residue_alternate_location.end(); it++)
    {
        residue_alternate_location_.push_back(*it);
    }
    selected_alternate_location_ = vector<bool>();
    for(vector<bool>::iterator it = selected_alternate_location.begin(); it != selected_alternate_location.end(); it++)
    {
        selected_alternate_location_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorAlternateResidue::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorAlternateResidue::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
string PdbPreprocessorAlternateResidue::GetResidueName()
{
    return residue_name_;
}
char PdbPreprocessorAlternateResidue::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
vector<char> PdbPreprocessorAlternateResidue::GetResidueAlternateLocation()
{
    return residue_alternate_location_;
}
vector<bool> PdbPreprocessorAlternateResidue::GetSelectedAlternateLocation()
{
    return selected_alternate_location_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorAlternateResidue::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorAlternateResidue::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorAlternateResidue::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorAlternateResidue::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}
void PdbPreprocessorAlternateResidue::SetResidueAlternateLocation(vector<char> residue_alternate_location)
{
    residue_alternate_location_.clear();
    for(vector<char>::iterator it = residue_alternate_location.begin(); it != residue_alternate_location.end(); it++)
    {
        residue_alternate_location_.push_back(*it);
    }
}
void PdbPreprocessorAlternateResidue::SetSelectedAlternateLocation(vector<bool> selected_alternate_location)
{
    selected_alternate_location_.clear();
    for(vector<bool>::iterator it = selected_alternate_location.begin(); it != selected_alternate_location.end(); it++)
    {
        selected_alternate_location_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorAlternateResidue::Print(ostream &out)
{
    out << "Residue name: " << residue_name_
         << ", Chain id: " << residue_chain_id_
         << ", Sequence number: " << residue_sequence_number_
         << ", Insertion code: " << residue_insertion_code_ << endl;
         for(unsigned int i = 0; i < residue_alternate_location_.size(); i++)
         {
             out << "\t Alternate Location: " << residue_alternate_location_[i] << " " << selected_alternate_location_[i] << endl;
         }

         out << endl;
}









