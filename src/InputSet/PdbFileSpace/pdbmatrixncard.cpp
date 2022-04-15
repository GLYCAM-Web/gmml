#include "../../../includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbMatrixNCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbMatrixNCard::PdbMatrixNCard() {}
PdbMatrixNCard::PdbMatrixNCard(std::string &line) {
    record_name_ = line.substr(0, 5);
    gmml::Trim(record_name_);
    if(line.substr(5, 1) == " ")
        n_ = gmml::iNotSet;
    else
        n_ = gmml::ConvertString<int>(line.substr(5, 1));
    if(line.substr(7, 3) == "   ")
        serial_number_ = gmml::iNotSet;
    else
        serial_number_ = gmml::ConvertString<int>(line.substr(7, 3));
    if(line.substr(10,10) == "          ")
        transfomration_vector_.SetX(gmml::dNotSet);
    else
        transfomration_vector_.SetX(gmml::ConvertString<double>(line.substr(10, 10)));
    if(line.substr(20, 10) == "          ")
        transfomration_vector_.SetY(gmml::dNotSet);
    else
        transfomration_vector_.SetY(gmml::ConvertString<double>(line.substr(20, 10)));
    if(line.substr(30, 10) == "          ")
        transfomration_vector_.SetZ(gmml::dNotSet);
    else
        transfomration_vector_.SetZ(gmml::ConvertString<double>(line.substr(30, 10)));
    if(line.substr(45, 10) == "          ")
        v_ = gmml::dNotSet;
    else
        v_ = gmml::ConvertString<double>(line.substr(45, 10));
    if(line.substr(59, 1) == " ")
        i_given_ = gmml::iNotSet;
    else
        i_given_ = gmml::ConvertString<int>(line.substr(59, 1));
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string PdbMatrixNCard::GetRecordName(){
    return record_name_;
}

int PdbMatrixNCard::GetN(){
    return n_;
}

int PdbMatrixNCard::GetSerialNumber(){
    return serial_number_;
}

GeometryTopology::Coordinate PdbMatrixNCard::GetTransformationVector(){
    return transfomration_vector_;
}

double PdbMatrixNCard::GetV(){
    return v_;
}

int PdbMatrixNCard::GetIGiven(){
    return i_given_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbMatrixNCard::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

void PdbMatrixNCard::SetN(int n){
    n_ = n;
}

void PdbMatrixNCard::SetSerialNumber(int serial_number){
    serial_number_ = serial_number;
}

void PdbMatrixNCard::SetTransformationVector(GeometryTopology::Coordinate transfomration_vector){
    transfomration_vector_ = transfomration_vector;
}

void PdbMatrixNCard::SetV(double v){
    v_ = v;
}

void PdbMatrixNCard::SetIGiven(int i_given){
    i_given_ = i_given;
}

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbMatrixNCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_;
    if(n_ != gmml::iNotSet)
        out << n_;
    else
        out << " ";
    out << ", Serial Number: ";
    if(serial_number_ != gmml::iNotSet)
        out << serial_number_;
    else
        out << " ";
    out << ", Transformation Vector: ";
    transfomration_vector_.Print(out);
    out << ", V: ";
    if(v_ != gmml::dNotSet)
        out << v_;
    else
        out << " ";
    out << ", I Given: ";
    if(i_given_ != gmml::iNotSet)
        out << i_given_;
    else
        out << " ";
    out << std::endl;
}
