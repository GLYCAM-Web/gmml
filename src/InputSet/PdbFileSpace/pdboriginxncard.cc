#include "../../../includes/InputSet/PdbFileSpace/pdboriginxncard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbOriginXnCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbOriginXnCard::PdbOriginXnCard():origin_() {}

PdbOriginXnCard::PdbOriginXnCard(std::stringstream& stream_block)
{
    std::string line;
    origin_=GeometryTopology::Coordinate();
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        record_name_ = line.substr(0,5);
        gmml::Trim(record_name_);
        if(line.substr(5,1) == " ")
            n_ = gmml::iNotSet;
        else
            n_ = gmml::ConvertString<int>(line.substr(5,1));
        if(line.substr(10,10) == "          ")
            origin_.SetX(gmml::dNotSet);
        else
            origin_.SetX( gmml::ConvertString<double>(line.substr(10,10)));
        if(line.substr(20, 10) == "          ")
            origin_.SetY(gmml::dNotSet);
        else
            origin_.SetY( gmml::ConvertString<double>(line.substr(20,10)));
        if(line.substr(30, 10) == "          ")
            origin_.SetZ(gmml::dNotSet);
        else
            origin_.SetZ( gmml::ConvertString<double>(line.substr(30,10)));
        if(line.substr(45, 10) == "          ")
            t_ = gmml::dNotSet;
        else
            t_ = gmml::ConvertString<double>(line.substr(45,10));

        getline(stream_block, line);
        temp = line;
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string PdbOriginXnCard::GetRecordName(){
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

void PdbOriginXnCard::SetRecordName(const std::string record_name){
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
void PdbOriginXnCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_;
    if(n_ != gmml::iNotSet)
        out << n_;
    else
        out << " ";
    out << ", Origin: ";
    origin_.Print(out);
    out << ", T: ";
    if(t_ != gmml::dNotSet)
        out << t_;
    else
        out << " ";
    out << std::endl;
}
