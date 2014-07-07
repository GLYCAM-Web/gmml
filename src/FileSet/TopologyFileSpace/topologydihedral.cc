
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyDihedral::TopologyDihedral() {}

TopologyDihedral::TopologyDihedral(vector<string> dihedral_atoms, vector<string> residue_names)
{
    dihedrals_.clear();
    for(vector<string>::iterator it = dihedral_atoms.begin(); it != dihedral_atoms.end(); it++)
    {
        dihedrals_.push_back(*it);
    }
    residue_names_.clear();
    for(vector<string>::iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
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
vector<string> TopologyDihedral::GetResidueNames()
{
    return residue_names_;
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
void TopologyDihedral::SetResidueNames(vector<string> residue_names)
{
    residue_names_.clear();
    for(vector<string>::iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyDihedral::Print(ostream &out)
{
    out << "Dihedral: " << residue_names_.at(0) << ":" << dihedrals_.at(0) << "-" << residue_names_.at(1) << ":" << dihedrals_.at(1) << "-" <<  residue_names_.at(3) << ":"  <<
           dihedrals_.at(2) << "-" << residue_names_.at(3) << ":" << dihedrals_.at(3) << endl;
    out << "\t ";
    dihedral_type_->Print(out);
    out << ", Including Hydrogen: ";
    if(including_hydrogen_)
        out << "YES";
    else
        out << "NO";
    out << ", Improper: ";
    if(is_improper_)
        out << "YES";
    else
        out << "NO";
    out << ", Ignored Group Interaction: ";
    if(ignored_group_interaction_)
        out << "YES";
    else
        out << "NO";
    out << endl;
}





