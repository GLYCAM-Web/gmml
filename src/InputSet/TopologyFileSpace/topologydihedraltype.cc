
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyDihedralType::TopologyDihedralType() {}

TopologyDihedralType::TopologyDihedralType(int index, double force_constant, double periodicity, double phase, double scee, double scnb) :
    index_(index), force_constant_(force_constant), periodicity_(periodicity), phase_(phase), scee_(scee), scnb_(scnb) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int TopologyDihedralType::GetIndex()
{
    return index_;
}
double TopologyDihedralType::GetPeriodicity()
{
    return periodicity_;
}
double TopologyDihedralType::GetPhase()
{
    return phase_;
}
double TopologyDihedralType::GetScee()
{
    return scee_;
}
double TopologyDihedralType::GetScnb()
{
    return scnb_;
}
double TopologyDihedralType::GetForceConstant()
{
    return force_constant_;
}


//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyDihedralType::SetIndex(int index)
{
    index_ = index;
}
void TopologyDihedralType::SetPeriodicity(double periodicity)
{
    periodicity_ = periodicity;
}
void TopologyDihedralType::SetPhase(double phase)
{
    phase_ = phase;
}
void TopologyDihedralType::SetScee(double scee)
{
    scee_ = scee;
}
void TopologyDihedralType::SetScnb(double scnb)
{
    scnb_ = scnb;
}
void TopologyDihedralType::SetForceConstant(double force_constant)
{
    force_constant_ = force_constant;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyDihedralType::Print(ostream &out)
{
    out << "Dihedral index: " << index_
         << ", Periodicity: " << periodicity_
         << ", Phase: " << phase_
         << ", SCEE: " << scee_
         << ", SCNB: " << scnb_
         << ", Force constant: " << force_constant_;
}



