#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp"
#include <iomanip>

using PdbPreprocessorSpace::PdbPreprocessorResidueInfo;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorResidueInfo::PdbPreprocessorResidueInfo() {}

PdbPreprocessorResidueInfo::PdbPreprocessorResidueInfo(
        std::string residue_name, char chain_id, int sequence_number,
        char residue_insertion_code, char residue_alternate_location, double residue_charge) :
    residue_chain_id_(chain_id), residue_sequence_number_(sequence_number), residue_name_(residue_name), 
    residue_insertion_code_(residue_insertion_code), residue_alternate_location_(residue_alternate_location),
    residue_charge_(residue_charge) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorResidueInfo::GetResidueChainId()
{
    return residue_chain_id_;
}
int PdbPreprocessorResidueInfo::GetResidueSequenceNumber()
{
    return residue_sequence_number_;
}
std::string PdbPreprocessorResidueInfo::GetResidueName()
{
    return residue_name_;
}
char PdbPreprocessorResidueInfo::GetResidueInsertionCode()
{
    return residue_insertion_code_;
}
char PdbPreprocessorResidueInfo::GetResidueAlternateLocation()
{
    return residue_alternate_location_;
}
double PdbPreprocessorResidueInfo::GetResidueCharge()
{
    return residue_charge_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorResidueInfo::SetResidueChainId(char residue_chain_id)
{
    residue_chain_id_ = residue_chain_id;
}
void PdbPreprocessorResidueInfo::SetResidueSequenceNumber(int residue_sequence_number)
{
    residue_sequence_number_ = residue_sequence_number;
}
void PdbPreprocessorResidueInfo::SetResidueName(const std::string residue_name)
{
    residue_name_ = residue_name;
}
void PdbPreprocessorResidueInfo::SetResidueInsertionCode(char residue_insertion_code)
{
    residue_insertion_code_ = residue_insertion_code;
}
void PdbPreprocessorResidueInfo::SetResidueAlternateLocation(char residue_alternate_location)
{
    residue_alternate_location_ = residue_alternate_location;
}
void PdbPreprocessorResidueInfo::SetResidueCharge(double residue_charge)
{
    residue_charge_ = residue_charge;
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorResidueInfo::Print(std::ostream &out)
{
    out << "Residue name: " << residue_name_
         << ", Chain id: " << residue_chain_id_
         << ", Sequence number: " << residue_sequence_number_
         << ", Insertion code: " << residue_insertion_code_
         << ", Alternate location: " << residue_alternate_location_
         << ", Residue charge: " << std::fixed << std::setprecision(2) << residue_charge_
         << std::endl;
}
