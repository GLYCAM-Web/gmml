
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"

using TopologyFileSpace::TopologyDihedral;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyDihedral::TopologyDihedral() {}

TopologyDihedral::TopologyDihedral(std::vector<std::string> dihedral_atoms, std::vector<std::string> residue_names)
{
    dihedrals_.clear();
    for(std::vector<std::string>::iterator it = dihedral_atoms.begin(); it != dihedral_atoms.end(); it++)
    {
        dihedrals_.push_back(*it);
    }
    residue_names_.clear();
    for(std::vector<std::string>::iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::vector<std::string> TopologyDihedral::GetDihedrals()
{
    return dihedrals_;
}
TopologyFileSpace::TopologyDihedralType* TopologyDihedral::GetDihedralType()
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
std::vector<std::string> TopologyDihedral::GetResidueNames()
{
    return residue_names_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyDihedral::SetDihedrals(std::vector<std::string> dihedrals)
{
    dihedrals_.clear();
    for(std::vector<std::string>::iterator it = dihedrals.begin(); it != dihedrals.end(); it++)
    {
        dihedrals_.push_back(*it);
    }
}
void TopologyDihedral::SetDihedralType(TopologyFileSpace::TopologyDihedralType* dihedral_type)
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
void TopologyDihedral::SetResidueNames(std::vector<std::string> residue_names)
{
    residue_names_.clear();
    for(std::vector<std::string>::iterator it = residue_names.begin(); it != residue_names.end(); it++)
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
void TopologyDihedral::Print(std::ostream &out)
{
    out << "Dihedral: " << residue_names_.at(0) << ":" << dihedrals_.at(0) << "-" << residue_names_.at(1) << ":" << dihedrals_.at(1) << "-" <<  residue_names_.at(3) << ":"  <<
           dihedrals_.at(2) << "-" << residue_names_.at(3) << ":" << dihedrals_.at(3) << std::endl;
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
    out << std::endl;
}
