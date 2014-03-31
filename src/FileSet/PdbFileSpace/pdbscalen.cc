#include "../../../includes/FileSet/PdbFileSpace/pdbscalen.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace Geometry;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbScaleN::PdbScaleN():scale_vector_() {}

PdbScaleN::PdbScaleN(stringstream& stream_block)
{
    string line;
    scale_vector_=Coordinate();
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        record_name_ = line.substr(0,5);
        Trim(record_name_);
        if(line.substr(5,1) == " ")
            n_ = iNotSet;
        else
            n_ = ConvertString<int>(line.substr(5,1));
        if(line.substr(10,10) == "          ")
            scale_vector_.SetX(dNotSet);
        else
            scale_vector_.SetX( ConvertString<double>(line.substr(10,10)));
        if(line.substr(20, 10) == "          ")
            scale_vector_.SetY(dNotSet);
        else
            scale_vector_.SetY( ConvertString<double>(line.substr(20,10)));
        if(line.substr(30, 10) == "          ")
            scale_vector_.SetZ(dNotSet);
        else
            scale_vector_.SetZ( ConvertString<double>(line.substr(30,10)));
        if(line.substr(45, 10) == "          ")
            u_ = dNotSet;
        else
            u_ = ConvertString<double>(line.substr(45,10));

        getline(stream_block, line);
        temp = line;
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
void PdbScaleN::Print(ostream &out)
{
    out << "Record Name: " << record_name_;
    if(n_ != iNotSet)
        out << n_;
    else
        out << " ";
    out << ", Origin: ";
    scale_vector_.Print(out);
    out << ", U: ";
    if(u_ != dNotSet)
        out << u_;
    else
        out << " ";
    out << endl;
}
