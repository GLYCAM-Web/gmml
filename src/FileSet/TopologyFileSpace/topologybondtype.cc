
#include "../../../includes/FileSet/TopologyFileSpace/topologybondtype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyBondType::TopologyBondType() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
vector<string> TopologyBondType::GetBondTypes()
{
    return bond_types_;
}
int TopologyBondType::GetIndex()
{
    return index_;
}
double TopologyBondType::GetForceConstant()
{
    return force_constant_;
}
double TopologyBondType::GetEquilibriumValue()
{
    return equilibrium_value_;
}


//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyBondType::SetBondTypes(vector<string> bond_types)
{   
    bond_types_.clear();
    for(vector<string>::iterator it = bond_types.begin(); it != bond_types.end(); it++)
    {
        bond_types_.push_back(*it);
    }
}
void TopologyBondType::SetIndex(int index)
{
    index_ = index;
}
void TopologyBondType::SetForceConstant(double force_constant)
{
    force_constant_ = force_constant;
}
void TopologyBondType::SetEquilibriumValue(double equilibrium_value)
{
    equilibrium_value_ = equilibrium_value;
}
//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyBondType::Print(ostream &out)
{
}


