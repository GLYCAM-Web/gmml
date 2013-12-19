#include "../../../includes/FileSet/PdbFileSpace/pdboriginxn.hpp"
//#include "../../../includes/Geometry/coordinate.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbOriginXn::PdbOriginXn() {}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbOriginXn::GetRecordName(){
    return record_name_;
}

int PdbOriginXn::GetN(){
    return n_;
}

Geometry::Coordinate PdbOriginXn::GetOrigin(){
    return origin_;
}

double PdbOriginXn::GetT(){
    return t_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbOriginXn::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbOriginXn::SetN(int n){
    n_ = n;
}

void PdbOriginXn::SetOrigin(Geometry::Coordinate origin){
    origin_ = origin;
}

void PdbOriginXn::SetT(double t){
    t_ = t;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////





