#include "../../../includes/InputSet/PdbFileSpace/pdbscalencard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbScaleNCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbScaleNCard::PdbScaleNCard():scale_vector_() {}

PdbScaleNCard::PdbScaleNCard(std::stringstream& stream_block)
{
    std::string line;
    scale_vector_=GeometryTopology::Coordinate();
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
            scale_vector_.SetX(gmml::dNotSet);
        else
            scale_vector_.SetX( gmml::ConvertString<double>(line.substr(10,10)));
        if(line.substr(20, 10) == "          ")
            scale_vector_.SetY(gmml::dNotSet);
        else
            scale_vector_.SetY( gmml::ConvertString<double>(line.substr(20,10)));
        if(line.substr(30, 10) == "          ")
            scale_vector_.SetZ(gmml::dNotSet);
        else
            scale_vector_.SetZ( gmml::ConvertString<double>(line.substr(30,10)));
        if(line.substr(45, 10) == "          ")
            u_ = gmml::dNotSet;
        else
            u_ = gmml::ConvertString<double>(line.substr(45,10));

        getline(stream_block, line);
        temp = line;
    }
}
//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string PdbScaleNCard::GetRecordName(){
    return record_name_;
}

int PdbScaleNCard::GetN(){
    return n_;
}

GeometryTopology::Coordinate PdbScaleNCard::GetScaleVector(){
    return scale_vector_;
}

double PdbScaleNCard::GetU(){
    return u_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbScaleNCard::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

void PdbScaleNCard::SetN(int n){
    n_ = n;
}

void PdbScaleNCard::SetScaleVector(GeometryTopology::Coordinate scale_vector){
    scale_vector_ = scale_vector;
}

void PdbScaleNCard::SetU(double u){
    u_ = u;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbScaleNCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_;
    if(n_ != gmml::iNotSet)
        out << n_;
    else
        out << " ";
    out << ", Origin: ";
    scale_vector_.Print(out);
    out << ", U: ";
    if(u_ != gmml::dNotSet)
        out << u_;
    else
        out << " ";
    out << std::endl;
}
