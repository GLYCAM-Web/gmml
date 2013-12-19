#include "../../../includes/FileSet/PdbFileSpace/pdbscalen.hpp"
//#include "../../../includes/Geometry/coordinate.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbScaleN::PdbScaleN() {}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbScaleN::GetRecordName(){
    return record_name_;
}

int PdbScaleN::GetN(){
    return n_;
}

Geometry::Coordinate PdbScaleN::GetScaleVector(){
    return scale_vector_;
}

double PdbScaleN::GetU(){
    return u_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbScaleN::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbScaleN::SetN(int n){
    n_ = n;
}

void PdbScaleN::SetScaleVector(Geometry::Coordinate scale_vector){
    scale_vector_ = scale_vector;
}

void PdbScaleN::SetU(double u){
    u_ = u;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////




