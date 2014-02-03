#include "../../../includes/FileSet/PdbFileSpace/pdbmatrixn.hpp"
//#include "../../../includes/Geometry/coordinate.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;
using namespace Geometry;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbMatrixN::PdbMatrixN() {}
PdbMatrixN::PdbMatrixN(string &line) {
    record_name_ = line.substr(0, 5);
    Trim(record_name_);
    n_ = ConvertString<int>(line.substr(5, 1));
    serial_number_ = ConvertString<int>(line.substr(7, 3));
    transfomration_vector_.SetX(ConvertString<double>(line.substr(10, 10)));
    transfomration_vector_.SetY(ConvertString<double>(line.substr(20, 10)));
    transfomration_vector_.SetZ(ConvertString<double>(line.substr(30, 10)));
    v_ = ConvertString<double>(line.substr(45, 10));
    i_given_ = ConvertString<int>(line.substr(59, 1));
}


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
void PdbMatrixN::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << n_ << ", Serial Number: " << serial_number_ << ", Transformation Vector: ";
    transfomration_vector_.Print(out);
    out << ", V: " << v_ << ", I Given: " << i_given_ << endl;
}
