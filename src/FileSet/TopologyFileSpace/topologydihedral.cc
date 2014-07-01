
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"

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
bool TopologyDihedral::GetIsImproper()
{
    return is_improper_;
}
bool TopologyDihedral::GetIgnoredGroupInteraction()
{
    return ignored_group_interaction_;
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
void TopologyDihedral::SetIsImproper(bool is_improper)
{
    is_improper_ = is_improper;
}
void TopologyDihedral::SetIgnoredGroupInteraction(bool ignored_group_interaction)
{
    ignored_group_interaction_ = ignored_group_interaction;
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





