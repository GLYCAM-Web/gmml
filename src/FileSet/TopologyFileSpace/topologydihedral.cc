
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyDihedral::TopologyDihedral() {}

TopologyDihedral::TopologyDihedral(vector<string> dihedral_atoms)
{
    dihedrals_.clear();
    for(vector<string>::iterator it = dihedral_atoms.begin(); it != dihedral_atoms.end(); it++)
    {
        dihedrals_.push_back(*it);
    }
}

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
bool TopologyDihedral::GetIncludingHydrogen()
{
    return including_hydrogen_;
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
void TopologyDihedral::SetIncludingHydrogen(bool including_hydrogen)
{
    including_hydrogen_ = including_hydrogen;
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





