
#include "../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"

using TopologyFileSpace::TopologyBondType;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyBondType::TopologyBondType() {}

TopologyBondType::TopologyBondType(int index, double force_constant, double equilibrium_value) :
    index_(index), force_constant_(force_constant), equilibrium_value_(equilibrium_value) {}

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
void TopologyBondType::Print(std::ostream &out)
{
    out << "Bond index: " << index_
         << ", Force constant: " << force_constant_
         << ", Equilibrium value: " << equilibrium_value_;
}
