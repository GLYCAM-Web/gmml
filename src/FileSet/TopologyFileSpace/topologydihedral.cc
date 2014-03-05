
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyDihedral::TopologyDihedral() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
vector<string> TopologyDihedral::GetDihedrals()
{
    return dihedrals_;
}
TopologyDihedralType* TopologyDihedral::GetDihedralType()
{
    return dihedral_type_;
}
bool TopologyDihedral::GetIsImproper()
{
    return is_improper_;
}
bool TopologyDihedral::GetIgnoredGroupInteraction()
{
    return ignored_group_interaction_;
}
TopologyPositionalDihedralType TopologyDihedral::GetPositionalType()
{
    return positional_type_;
}
TopologyRotationalDihedralType TopologyDihedral::GetRotationalType()
{
    return rotational_type_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyDihedral::SetDihedrals(vector<string> dihedrals)
{
    dihedrals_.clear();
    for(vector<string>::iterator it = dihedrals.begin(); it != dihedrals.end(); it++)
    {
        dihedrals_.push_back(*it);
    }
}
void TopologyDihedral::SetDihedralType(TopologyDihedralType* dihedral_type)
{
    dihedral_type_ = dihedral_type;
}
void TopologyDihedral::SetIsImproper(bool is_improper)
{
    is_improper_ = is_improper;
}
void TopologyDihedral::SetIgnoredGroupInteraction(bool ignored_group_interaction)
{
    ignored_group_interaction_ = ignored_group_interaction;
}
void TopologyDihedral::SetPositionalType(TopologyPositionalDihedralType positional_type)
{
    positional_type_ = positional_type;
}
void TopologyDihedral::SetRotationalType(TopologyRotationalDihedralType rotational_type)
{
    rotational_type_ = rotational_type;
}
//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyDihedral::Print(ostream &out)
{
}





