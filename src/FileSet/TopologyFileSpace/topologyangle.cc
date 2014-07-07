
#include "../../../includes/FileSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyangletype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAngle::TopologyAngle() {}

TopologyAngle::TopologyAngle(vector<string> angle_atoms, vector<string> residue_names)
{
    angles_.clear();
    for(vector<string>::iterator it = angle_atoms.begin(); it != angle_atoms.end(); it++)
    {
        angles_.push_back(*it);
    }
    residue_names_.clear();
    for(vector<string>::iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
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
vector<string> TopologyAngle::GetResidueNames()
{
    return residue_names_;
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
void TopologyAngle::SetResidueNames(vector<string> residue_names)
{
    residue_names_.clear();
    for(vector<string>::iterator it = residue_names.begin(); it != residue_names.end(); it++)
    {
        residue_names_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAngle::Print(ostream &out)
{
    out << "Angle: " << residue_names_.at(0) << ":" << angles_.at(0) << "-" << residue_names_.at(1) << ":" << angles_.at(1) << "-" << residue_names_.at(2) << ":" << angles_.at(2) << endl;
    out << "\t ";
    angle_type_->Print(out);
    out << ", Including Hydrogen: ";
    if(including_hydrogen_)
        out << "YES";
    else
        out << "NO";
    out << endl;
}




