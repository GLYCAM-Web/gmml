
#include "../../../includes/FileSet/TopologyFileSpace/topologyangletype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAngleType::TopologyAngleType() {}

TopologyAngleType::TopologyAngleType(int index, double force_constant, double equilibrium_value) :
    index_(index), force_constant_(force_constant), equilibrium_value_(equilibrium_value) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int TopologyAngleType::GetIndex()
{
    return index_;
}
double TopologyAngleType::GetForceConstant()
{
    return force_constant_;
}
double TopologyAngleType::GetEquilibriumValue()
{
    return equilibrium_value_;
}


//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAngleType::SetIndex(int index)
{
    index_ = index;
}
void TopologyAngleType::SetForceConstant(double force_constant)
{
    force_constant_ = force_constant;
}
void TopologyAngleType::SetEquilibriumValue(double equilibrium_value)
{
    equilibrium_value_ = equilibrium_value;
}
//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAngleType::Print(ostream &out)
{
}



