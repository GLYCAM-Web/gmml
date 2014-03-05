
#include "../../../includes/FileSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyResidue::TopologyResidue() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string TopologyResidue::GetResidueName()
{
    return residue_name_;
}
TopologyResidue::TopologyAtomMap TopologyResidue::GetAtoms()
{
    return atoms_;
}
TopologyResidue::TopologyBondMap TopologyResidue::GetBonds()
{
    return bonds_;
}
TopologyResidue::TopologyAngleMap TopologyResidue::GetAngles()
{
    return angles_;
}
TopologyResidue::TopologyDihedralMap TopologyResidue::GetDihedrals()
{
    return dihedrals_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyResidue::SetResidueName(const string residue_name)
{
    residue_name_ = residue_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyResidue::Print(ostream &out)
{
}





