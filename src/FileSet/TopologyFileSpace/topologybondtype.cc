
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


