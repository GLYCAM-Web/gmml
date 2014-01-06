#include "../../../includes/FileSet/PdbFileSpace/pdbscalen.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace Geometry;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbScaleN::PdbScaleN():scale_vector_() {}

PdbScaleN::PdbScaleN(istringstream& stream_block)
{
    string line;
    scale_vector_=Coordinate();
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        record_name_ = line.substr(0,6);
        n_ = ConvertString<int>(line.substr(5,1));
        scale_vector_.SetX( ConvertString<double>(line.substr(10,10)));
        scale_vector_.SetY( ConvertString<double>(line.substr(20,10)));
        scale_vector_.SetZ( ConvertString<double>(line.substr(30,10)));
        u_ = ConvertString<double>(line.substr(45,10));

        getline(stream_block, line);
    }
}
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




