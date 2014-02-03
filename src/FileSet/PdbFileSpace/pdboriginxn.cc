#include "../../../includes/FileSet/PdbFileSpace/pdboriginxn.hpp"
//#include "../../../includes/Geometry/coordinate.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace Geometry;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbOriginXn::PdbOriginXn():origin_() {}

PdbOriginXn::PdbOriginXn(stringstream& stream_block)
{
    string line;
    origin_=Coordinate();
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        record_name_ = line.substr(0,5);
        Trim(record_name_);
        n_ = ConvertString<int>(line.substr(5,1));
        origin_.SetX( ConvertString<double>(line.substr(10,10)));
        origin_.SetY( ConvertString<double>(line.substr(20,10)));
        origin_.SetZ( ConvertString<double>(line.substr(30,10)));
        t_ = ConvertString<double>(line.substr(45,10));

        getline(stream_block, line);
        temp = line;
    }
}


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
void PdbOriginXn::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << n_ << ", Origin: ";
    origin_.Print(out);
    out << ", T: " << t_ << endl;
}
