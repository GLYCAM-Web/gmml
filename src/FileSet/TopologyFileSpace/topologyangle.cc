
#include "../../../includes/FileSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangletype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAngle::TopologyAngle() {}

TopologyAngle::TopologyAngle(vector<string> angle_atoms)
{
    angles_.clear();
    for(vector<string>::iterator it = angle_atoms.begin(); it != angle_atoms.end(); it++)
    {
        angles_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
vector<string> TopologyAngle::GetAngles()
{
    return angles_;
}
TopologyAngleType* TopologyAngle::GetAngleType()
{
    return angle_type_;
}
bool TopologyAngle::GetIncludingHydrogen()
{
    return including_hydrogen_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAngle::SetAngles(vector<string> angles)
{
    angles_.clear();
    for(vector<string>::iterator it = angles.begin(); it != angles.end(); it++)
    {
        angles_.push_back(*it);
    }
}
void TopologyAngle::SetAnlgeType(TopologyAngleType* angle_type)
{
    angle_type_ = angle_type;
}
void TopologyAngle::SetIncludingHydrogen(bool including_hydrogen)
{
    including_hydrogen_ = including_hydrogen;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAngle::Print(ostream &out)
{
    out << "Angle: " << angles_.at(0) << "-" << angles_.at(1) << "-" << angles_.at(2) << endl;
    out << "\t ";
    angle_type_->Print(out);
    out << "\t Including Hydrogen: ";
    if(including_hydrogen_)
        out << "YES";
    else
        out << "NO";
    out << endl;
}




