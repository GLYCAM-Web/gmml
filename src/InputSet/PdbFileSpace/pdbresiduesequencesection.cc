// Author: Alireza Khatamian

#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequencesection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbResidueSequenceSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbResidueSequenceSection::PdbResidueSequenceSection() : record_name_("SEQRES") {}

PdbResidueSequenceSection::PdbResidueSequenceSection(const std::string &record_name) : record_name_(record_name) {}

PdbResidueSequenceSection::PdbResidueSequenceSection(std::stringstream &stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
        std::stringstream residue_sequence_block;
        residue_sequence_block << line << std::endl;
        char chain_id = gmml::ConvertString<char>(line.substr(11,1));

        getline(stream_block, line);
        temp = line;

        while (!gmml::Trim(temp).empty() && gmml::ConvertString<char>(line.substr(11,1)) == chain_id){
            residue_sequence_block << line << std::endl;
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
std::string PdbResidueSequenceSection::GetRecordName()
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
void PdbResidueSequenceSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbResidueSequenceSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "============= Residue Sequence Chain =============" << std::endl;
    for(PdbResidueSequenceSection::ResidueSequenceCardMap::iterator it = residue_sequence_chains_.begin(); it != residue_sequence_chains_.end(); it++)
    {
        out << "Chain ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
