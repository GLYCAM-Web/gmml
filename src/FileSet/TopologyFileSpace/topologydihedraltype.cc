
#include "../../../includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyDihedralType::TopologyDihedralType() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
vector<string> TopologyDihedralType::GetDihedralTypes()
{
    return dihedral_types_;
}
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


//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyDihedralType::SetDihedralType(vector<string> dihedral_types)
{
    dihedral_types_.clear();
    for(vector<string>::iterator it = dihedral_types.begin(); it != dihedral_types.end(); it++)
    {
        dihedral_types_.push_back(*it);
    }
}
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
//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyDihedralType::Print(ostream &out)
{
}



