#include "../../includes/MolecularModeling/oligosaccharidedetectionatom.hpp"

using MolecularModeling::OligoSaccharideDetectionAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
OligoSaccharideDetectionAtom::OligoSaccharideDetectionAtom() : IsCycle_(false), IsSideChain_(false), IsAnomericCarbon_(false) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
bool OligoSaccharideDetectionAtom::GetIsCycle()
{
    return IsCycle_;
}
bool OligoSaccharideDetectionAtom::GetIsSideChain()
{
    return IsSideChain_;
}
bool OligoSaccharideDetectionAtom::GetIsAnomericCarbon()
{
    return IsAnomericCarbon_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void OligoSaccharideDetectionAtom::SetIsCycle(bool is_cycle)
{
    IsCycle_ = is_cycle;
}
void OligoSaccharideDetectionAtom::SetIsSideChain(bool is_side_chain)
{
    IsSideChain_ = is_side_chain;
}
void OligoSaccharideDetectionAtom::SetIsAnomericCarbon(bool is_anomeric_carbon)
{
    IsAnomericCarbon_ = is_anomeric_carbon;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void OligoSaccharideDetectionAtom::Print(std::ostream &out)
{
	out << "";
}
