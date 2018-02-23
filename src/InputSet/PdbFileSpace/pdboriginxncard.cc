#include "../../../includes/InputSet/PdbFileSpace/pdboriginxncard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace GeometryTopology;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbOriginXnCard::PdbOriginXnCard():origin_() {}

PdbOriginXnCard::PdbOriginXnCard(stringstream& stream_block)
{
    string line;
    origin_=Coordinate();
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
            origin_.SetX(dNotSet);
        else
            origin_.SetX( ConvertString<double>(line.substr(10,10)));
        if(line.substr(20, 10) == "          ")
            origin_.SetY(dNotSet);
        else
            origin_.SetY( ConvertString<double>(line.substr(20,10)));
        if(line.substr(30, 10) == "          ")
            origin_.SetZ(dNotSet);
        else
            origin_.SetZ( ConvertString<double>(line.substr(30,10)));
        if(line.substr(45, 10) == "          ")
            t_ = dNotSet;
        else
            t_ = ConvertString<double>(line.substr(45,10));

        getline(stream_block, line);
        temp = line;
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbOriginXnCard::GetRecordName(){
    return record_name_;
}

int PdbOriginXnCard::GetN(){
    return n_;
}

GeometryTopology::Coordinate PdbOriginXnCard::GetOrigin(){
    return origin_;
}

double PdbOriginXnCard::GetT(){
    return t_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbOriginXnCard::SetRecordName(const string record_name){
    record_name_ = record_name;
}

void PdbOriginXnCard::SetN(int n){
    n_ = n;
}

void PdbOriginXnCard::SetOrigin(GeometryTopology::Coordinate origin){
    origin_ = origin;
}

void PdbOriginXnCard::SetT(double t){
    t_ = t;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbOriginXnCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_;
    if(n_ != iNotSet)
        out << n_;
    else
        out << " ";
    out << ", Origin: ";
    origin_.Print(out);
    out << ", T: ";
    if(t_ != dNotSet)
        out << t_;
    else
        out << " ";
    out << endl;
}
