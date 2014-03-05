
#include "../../../includes/FileSet/TopologyFileSpace/topologyangletype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAngleType::TopologyAngleType() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
vector<string> TopologyAngleType::GetAngleTypes()
{
    return angle_types_;
}
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
void TopologyAngleType::SetAngleTypes(vector<string> angle_types)
{   
    angle_types_.clear();
    for(vector<string>::iterator it = angle_types.begin(); it != angle_types.end(); it++)
    {
        angle_types_.push_back(*it);
    }
}
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



