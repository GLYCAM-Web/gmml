#include "../../includes/MolecularModeling/oligosaccharidedetectionatom.hpp"

using MolecularModeling::OligoSaccharideDetectionAtom;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
MolecularModeling::OligoSaccharideDetectionAtom::OligoSaccharideDetectionAtom() : IsCycle_(false), IsSideChain_(false), IsAnomericCarbon_(false), naming_("") {}

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
std::string MolecularModeling::OligoSaccharideDetectionAtom::GetNaming()
{
    return naming_;
}
OligoSaccharideDetectionAtom::OligosaccharidePropertyTags MolecularModeling::OligoSaccharideDetectionAtom::GetAllOligosaccharidePropertyTags()
{
    return oligosaccharide_properties_;
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
void MolecularModeling::OligoSaccharideDetectionAtom::SetNaming(std::string naming)
{
    naming_ = naming;
}
void MolecularModeling::OligoSaccharideDetectionAtom::AddOligosaccharidePropertyTag (std::pair<std::string, std::string> role_tag_pair)
{
    oligosaccharide_properties_.insert(role_tag_pair);
}
void MolecularModeling::OligoSaccharideDetectionAtom::RemoveOligosaccharidePropertyTag (std::pair<std::string, std::string> role_tag_pair)
{
    std::string role = role_tag_pair.first;
    std::string tag = role_tag_pair.second;
    for (OligosaccharidePropertyTags::iterator it = oligosaccharide_properties_.begin(); it != oligosaccharide_properties_.end(); it++ ){
	std::string map_role = it->first;
	std::string map_tag = it->second;
	if (map_role == role && map_tag == tag){
	    oligosaccharide_properties_.erase(it);
	}
    }
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void MolecularModeling::OligoSaccharideDetectionAtom::Print(std::ostream &out)
{
	out << "";
}
