#include "includes/InputSet/PdbFile/pdbResidue.hpp"

using pdb::Residue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Residue::Residue() {}

Residue::Residue(std::vector<AtomRecord*> atomRecords){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
//char Residue::GetResidueChainId()
//{
//    return residue_chain_id_;
//}
//std::string Residue::GetResidueName()
//{
//    return residue_name_;
//}
//int Residue::GetResidueSequenceNumber()
//{
//    return residue_sequence_number_;
//}
//char Residue::GetResidueInsertionCode()
//{
//    return residue_insertion_code_;
//}
//char Residue::GetResidueAlternateLocation()
//{
//    return residue_alternate_location_;
//}
//
////////////////////////////////////////////////////////////
////                          MUTATOR                     //
////////////////////////////////////////////////////////////
//void Residue::SetResidueChainId(char residue_chain_id)
//{
//    residue_chain_id_ = residue_chain_id;
//}
//void Residue::SetResidueName(const std::string residue_name)
//{
//    residue_name_ = residue_name;
//}
//void Residue::SetResidueSequenceNumber(int residue_sequence_number)
//{
//    residue_sequence_number_ = residue_sequence_number;
//}
//void Residue::SetResidueInsertionCode(char residue_insertion_code)
//{
//    residue_insertion_code_ = residue_insertion_code;
//}
//void Residue::SetResidueAlternateLocation(char residue_alternate_location)
//{
//    residue_alternate_location_ = residue_alternate_location;
//}
//
////////////////////////////////////////////////////////////
////                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////////
//void Residue::Print(std::ostream &out) const
//{
//    out << "Residue name: " << residue_name_
//         << ", Chain id: " << residue_chain_id_
//         << ", Sequence number: " << residue_sequence_number_
//         << ", Insertion code: " << residue_insertion_code_
//         << ", Alternate location: " << residue_alternate_location_ << std::endl;
//}
