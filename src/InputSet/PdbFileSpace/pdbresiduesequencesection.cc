// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequencesection.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueSequenceSection::PdbResidueSequenceSection() : record_name_("SEQRES") {}

PdbResidueSequenceSection::PdbResidueSequenceSection(const string &record_name) : record_name_(record_name) {}

PdbResidueSequenceSection::PdbResidueSequenceSection(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
        stringstream residue_sequence_block;
        residue_sequence_block << line << endl;
        char chain_id = ConvertString<char>(line.substr(11,1));

        getline(stream_block, line);
        temp = line;

        while (!Trim(temp).empty() && ConvertString<char>(line.substr(11,1)) == chain_id){
            residue_sequence_block << line << endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbResidueSequenceCard* residue_sequence = new PdbResidueSequenceCard(residue_sequence_block);
        residue_sequence_chains_[chain_id] = residue_sequence;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbResidueSequenceSection::GetRecordName()
{
    return record_name_;
}

PdbResidueSequenceSection::ResidueSequenceCardMap PdbResidueSequenceSection::GetResidueSequenceChain()
{
    return residue_sequence_chains_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueSequenceSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueSequenceSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "============= Residue Sequence Chain =============" << endl;
    for(PdbResidueSequenceSection::ResidueSequenceCardMap::iterator it = residue_sequence_chains_.begin(); it != residue_sequence_chains_.end(); it++)
    {
        out << "Chain ID: " << (it)->first << endl;
        (it)->second->Print();
        out << endl;
    }
}
