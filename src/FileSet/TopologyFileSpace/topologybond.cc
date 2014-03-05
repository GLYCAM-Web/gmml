
#include "../../../includes/FileSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologybondtype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyBond::TopologyBond() {}

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
TopologyPositionalBondType TopologyBond::GetPositionalType()
{
    return positional_type_;
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
void TopologyBond::SetPositionalType(TopologyPositionalBondType positional_type)
{
    positional_type_ = positional_type;
}
//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyBond::Print(ostream &out)
{
}



