#include "../../../includes/MolecularModeling/residueproperties.hpp"

using MolecularModeling::ResidueProperties;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

ResidueProperties::ResidueProperties() :  is_residue_solvent_(false){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////


bool ResidueProperties::GetIsResidueSolvent()
{
    return is_residue_solvent_;
}


//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

void ResidueProperties::SetIsResidueSolvent(bool is_residue_solvent)
{
    is_residue_solvent_ = is_residue_solvent;
}


//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////


void ResidueProperties::Print(std::ostream &out)
{
     out << "------------------------ Residue Properties --------------------------" << std::endl;
     out << "Is Residue Solvent : " << is_residue_solvent_ << std::endl;
}
