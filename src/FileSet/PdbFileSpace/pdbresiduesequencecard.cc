// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbresiduesequencecard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbresiduesequence.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueSequenceCard::PdbResidueSequenceCard() : record_name_("SEQRES") {}

PdbResidueSequenceCard::PdbResidueSequenceCard(const string &record_name) : record_name_(record_name) {}

PdbResidueSequenceCard::PdbResidueSequenceCard(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }
        stringstream residue_sequence_block;
        residue_sequence_block << line << endl;
        char chain_id = ConvertString<char>(line.substr(11,1));

        getline(stream_block, line);

        while (!Trim(line).empty() && ConvertString<char>(line.substr(11,1)) == chain_id){
            residue_sequence_block << line << endl;
            getline(stream_block, line);
        }
        PdbResidueSequence* residue_sequence = new PdbResidueSequence(residue_sequence_block);
        residue_sequence_chains_[chain_id] = residue_sequence;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbResidueSequenceCard::GetRecordName()
{
    return record_name_;
}

PdbResidueSequenceCard::ResidueSequenceMap PdbResidueSequenceCard::GetResidueSequenceChain()
{
    return residue_sequence_chains_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbResidueSequenceCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

