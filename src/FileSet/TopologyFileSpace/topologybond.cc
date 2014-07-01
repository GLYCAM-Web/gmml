
#include "../../../includes/FileSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybondtype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyBond::TopologyBond() {}

TopologyBond::TopologyBond(vector<string> bonds)
{
    bonds_.clear();
    for(vector<string>::iterator it = bonds.begin(); it != bonds.end(); it++)
    {
        bonds_.push_back(*it);
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

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyBond::Print(ostream &out)
{
    out << "Bond: " << bonds_.at(0) << "-" << bonds_.at(1) << endl;
    out << "\t ";
    bond_type_->Print(out);
    out << "\t Including Hydrogen: ";
    if(including_hydrogen_)
        out << "YES";
    else
        out << "NO";
    out << endl;
}



