#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp"

using PdbPreprocessorSpace::PdbPreprocessorDisulfideBond;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessorDisulfideBond::PdbPreprocessorDisulfideBond() {}

PdbPreprocessorDisulfideBond::PdbPreprocessorDisulfideBond(
        char residue_chain_id_1, char residue_chain_id_2, int residue_sequence_number_1,
        int residue_sequence_number_2, double distance, bool is_bonded,
        char residue_insertion_code_1, char residue_insertion_code_2,
        char residue_alternate_location_1, char residue_alternate_location_2,
        int sulfur_atom_serial_number_1, int sulfur_atom_serial_number_2) :
    residue_chain_id_1_(residue_chain_id_1), residue_chain_id_2_(residue_chain_id_2), residue_sequence_number_1_(residue_sequence_number_1),
    residue_sequence_number_2_(residue_sequence_number_2), distance_(distance), is_bonded_(is_bonded), residue_insertion_code_1_(residue_insertion_code_1),
    residue_insertion_code_2_(residue_insertion_code_2), residue_alternate_location_1_(residue_alternate_location_1),
    residue_alternate_location_2_(residue_alternate_location_2), sulfur_atom_serial_number_1_(sulfur_atom_serial_number_1),
    sulfur_atom_serial_number_2_(sulfur_atom_serial_number_2){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
char PdbPreprocessorDisulfideBond::GetResidueChainId1()
{
    return residue_chain_id_1_;
}
char PdbPreprocessorDisulfideBond::GetResidueChainId2()
{
    return residue_chain_id_2_;
}
int PdbPreprocessorDisulfideBond::GetResidueSequenceNumber1()
{
    return residue_sequence_number_1_;
}
int PdbPreprocessorDisulfideBond::GetResidueSequenceNumber2()
{
    return residue_sequence_number_2_;
}
double PdbPreprocessorDisulfideBond::GetDistance()
{
    return distance_;
}
bool PdbPreprocessorDisulfideBond::GetIsBonded()
{
    return is_bonded_;
}
char PdbPreprocessorDisulfideBond::GetResidueInsertionCode1()
{
    return residue_insertion_code_1_;
}
char PdbPreprocessorDisulfideBond::GetResidueInsertionCode2()
{
    return residue_insertion_code_2_;
}
char PdbPreprocessorDisulfideBond::GetResidueAlternateLocation1()
{
    return residue_alternate_location_1_;
}
char PdbPreprocessorDisulfideBond::GetResidueAlternateLocation2()
{
    return residue_alternate_location_2_;
}
int PdbPreprocessorDisulfideBond::GetSulfurAtomSerialNumber1()
{
    return sulfur_atom_serial_number_1_;
}
int PdbPreprocessorDisulfideBond::GetSulfurAtomSerialNumber2()
{
    return sulfur_atom_serial_number_2_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessorDisulfideBond::SetResidueChainId1(char residue_chain_id_1)
{
    residue_chain_id_1_ = residue_chain_id_1;
}
void PdbPreprocessorDisulfideBond::SetResidueChainId2(char residue_chain_id_2)
{
    residue_chain_id_2_ = residue_chain_id_2;
}
void PdbPreprocessorDisulfideBond::SetResidueSequenceNumber1(int residue_sequence_number_1)
{
    residue_sequence_number_1_ = residue_sequence_number_1;
}
void PdbPreprocessorDisulfideBond::SetResidueSequenceNumber2(int residue_sequence_number_2)
{
    residue_sequence_number_2_ = residue_sequence_number_2;
}
void PdbPreprocessorDisulfideBond::SetDistance(double distance)
{
    distance_ = distance;
}
void PdbPreprocessorDisulfideBond::SetIsBonded(bool is_bonded)
{
    is_bonded_ = is_bonded;
}
void PdbPreprocessorDisulfideBond::SetResidueInsertionCode1(char residue_insertion_code_1)
{
    residue_insertion_code_1_ = residue_insertion_code_1;
}
void PdbPreprocessorDisulfideBond::SetResidueInsertionCode2(char residue_insertion_code_2)
{
    residue_insertion_code_2_ = residue_insertion_code_2;
}
void PdbPreprocessorDisulfideBond::SetResidueAlternateLocation1(char residue_alternate_location_1)
{
    residue_alternate_location_1_ = residue_alternate_location_1;
}
void PdbPreprocessorDisulfideBond::SetResidueAlternateLocation2(char residue_alternate_location_2)
{
    residue_alternate_location_2_ = residue_alternate_location_2;
}
void PdbPreprocessorDisulfideBond::SetSulfurAtomSerialNumber1(int sulfur_atom_serial_number_1)
{
    sulfur_atom_serial_number_1_ = sulfur_atom_serial_number_1;
}
void PdbPreprocessorDisulfideBond::SetSulfurAtomSerialNumber2(int sulfur_atom_serial_number_2)
{
    sulfur_atom_serial_number_2_ = sulfur_atom_serial_number_2;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessorDisulfideBond::Print(std::ostream &out)
{
    out << "[Chain id: " << residue_chain_id_1_
         << ", Sequence number: " << residue_sequence_number_1_
         << ", Insertion code: " << residue_insertion_code_1_
         << ", Alternate location: " << residue_alternate_location_1_
         << "]--[Chain id: " << residue_chain_id_2_
         << ", Sequence number: " << residue_sequence_number_2_
         << ", Insertion code: " << residue_insertion_code_2_
         << ", Alternate location: " << residue_alternate_location_2_
         << "] => Distance: " << distance_
         << std::endl;
}
