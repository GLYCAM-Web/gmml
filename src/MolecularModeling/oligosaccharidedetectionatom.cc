#include "../../includes/MolecularModeling/oligosaccharidedetectionatom.hpp"

using MolecularModeling::OligoSaccharideDetectionAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
MolecularModeling::OligoSaccharideDetectionAtom::OligoSaccharideDetectionAtom() : IsCycle_(false), IsSideChain_(false), IsAnomericCarbon_(false) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
bool MolecularModeling::OligoSaccharideDetectionAtom::GetIsCycle()
{
    return IsCycle_;
}
bool MolecularModeling::OligoSaccharideDetectionAtom::GetIsSideChain()
{
    return IsSideChain_;
}
bool MolecularModeling::OligoSaccharideDetectionAtom::GetIsAnomericCarbon()
{
    return IsAnomericCarbon_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void MolecularModeling::OligoSaccharideDetectionAtom::SetIsCycle(bool is_cycle)
{
    IsCycle_ = is_cycle;
}
void MolecularModeling::OligoSaccharideDetectionAtom::SetIsSideChain(bool is_side_chain)
{
    IsSideChain_ = is_side_chain;
}
void MolecularModeling::OligoSaccharideDetectionAtom::SetIsAnomericCarbon(bool is_anomeric_carbon)
{
    IsAnomericCarbon_ = is_anomeric_carbon;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void MolecularModeling::OligoSaccharideDetectionAtom::Print(std::ostream &out)
{
	out << "";
}
