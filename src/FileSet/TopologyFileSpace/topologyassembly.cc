
#include "../../../includes/FileSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAssembly::TopologyAssembly() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string TopologyAssembly::GetAssemblyName()
{
    return assembly_name_;
}
TopologyAssembly::TopologyResidueMap TopologyAssembly::GetResidues()
{
    return residues_;
}
TopologyAssembly::TopologyBondMap TopologyAssembly::GetBonds()
{
    return bonds_;
}
TopologyAssembly::TopologyAngleMap TopologyAssembly::GetAngles()
{
    return angles_;
}
TopologyAssembly::TopologyDihedralMap TopologyAssembly::GetDihedrals()
{
    return dihedrals_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAssembly::SetAssemblyName(const string assembly_name)
{
    assembly_name_ = assembly_name;
}
void TopologyAssembly::SetResidues(TopologyResidueMap residues)
{
    residues_.clear();
    for(TopologyResidueMap::iterator it = residues.begin(); it != residues.end(); it++)
    {
        TopologyResidue* residue = (*it).second;
        residues_[residue->GetResidueName()] = residue;
    }
}
void TopologyAssembly::SetBonds(TopologyBondMap bonds)
{
    bonds_.clear();
    for(TopologyBondMap::iterator it = bonds.begin(); it != bonds.end(); it++)
    {
        TopologyBond* bond = (*it).second;
        bonds_[bond->GetBonds()] = bond;
    }
}
void TopologyAssembly::SetAngles(TopologyAngleMap angles)
{
    angles_.clear();
    for(TopologyAngleMap::iterator it = angles.begin(); it != angles.end(); it++)
    {
        TopologyAngle* angle = (*it).second;
        angles_[angle->GetAngles()] = angle;
    }
}
void TopologyAssembly::SetDihedrals(TopologyDihedralMap dihedrals)
{
    dihedrals_.clear();
    for(TopologyDihedralMap::iterator it = dihedrals.begin(); it != dihedrals.end(); it++)
    {
        TopologyDihedral* dihedral = (*it).second;
        dihedrals_[dihedral->GetDihedrals()] = dihedral;
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAssembly::Print(ostream &out)
{
}

