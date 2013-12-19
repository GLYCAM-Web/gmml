#include "../../../includes/FileSet/PdbFileSpace/pdbmatrixn.hpp"
//#include "../../../includes/Geometry/coordinate.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbMatrixN::PdbMatrixN() {}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbMatrixN::GetRecordName(){
    return record_name_;
}

int PdbMatrixN::GetN(){
    return n_;
}

int PdbMatrixN::GetSerialNumber(){
    return serial_number_;
}

Geometry::Coordinate PdbMatrixN::GetTransformationVector(){
    return transfomration_vector_;
}

double PdbMatrixN::GetV(){
    return v_;
}

int PdbMatrixN::GetIGiven(){
    return i_given_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbMatrixN::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbMatrixN::SetN(int n){
    n_ = n;
}

void PdbMatrixN::SetSerialNumber(int serial_number){
    serial_number_ = serial_number;
}

void PdbMatrixN::SetTransformationVector(Geometry::Coordinate transfomration_vector){
    transfomration_vector_ = transfomration_vector;
}

void PdbMatrixN::SetV(double v){
    v_ = v;
}

void PdbMatrixN::SetIGiven(int i_given){
    i_given_ = i_given;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////



