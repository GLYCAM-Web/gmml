
#include "../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyBond::TopologyBond() {}

TopologyBond::TopologyBond(vector<string> bonds, vector<string> residue_names)
{
    bonds_.clear();
    for(vector<string>::iterator it = bonds.begin(); it != bonds.end(); it++)
    {
        bonds_.push_back(*it);
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
vector<string> TopologyBond::GetBonds()
{
    return bonds_;
}
TopologyBondType* TopologyBond::GetBondType()
{
    return bond_type_;
}
bool TopologyBond::GetIncludingHydrogen()
{
    return including_hydrogen_;
}
vector<string> TopologyBond::GetResidueNames()
{
    return residue_names_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyBond::SetBonds(vector<string> bonds)
{
    bonds_.clear();
    for(vector<string>::iterator it = bonds.begin(); it != bonds.end(); it++)
    {
        bonds_.push_back(*it);
    }
}
void TopologyBond::SetBondType(TopologyBondType* bond_type)
{
    bond_type_ = bond_type;
}
void TopologyBond::SetIncludingHydrogen(bool including_hydrogen)
{
    including_hydrogen_ = including_hydrogen;
}
void TopologyBond::SetResidueNames(vector<string> residue_names)
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
void TopologyBond::Print(ostream &out)
{
    out << "Bond: " << residue_names_.at(0) << ":" << bonds_.at(0) << "-" << residue_names_.at(1) << ":" << bonds_.at(1) << endl;
    out << "\t ";
    bond_type_->Print(out);
    out << ", Including Hydrogen: ";
    if(including_hydrogen_)
        out << "YES";
    else
        out << "NO";
    out << endl;
}



