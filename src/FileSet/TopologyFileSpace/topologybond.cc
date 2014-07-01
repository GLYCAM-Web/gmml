
#include "../../../includes/FileSet/TopologyFileSpace/topologybond.hpp"

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
}



