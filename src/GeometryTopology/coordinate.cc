#include <math.h>

#include "../../includes/GeometryTopology/coordinate.hpp"
#include "../../includes/common.hpp"
#include "../../includes/utils.hpp"

using namespace GeometryTopology;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
Coordinate::Coordinate() : x_(0.0), y_(0.0), z_(0.0) {}

Coordinate::Coordinate(double x, double y, double z) : x_(x), y_(y), z_(z) {}

Coordinate::Coordinate(const Coordinate &coordinate) : x_(coordinate.x_), y_(coordinate.y_), z_(coordinate.z_) {}

Coordinate::Coordinate(Coordinate* coordinate) : x_(coordinate->x_), y_(coordinate->y_), z_(coordinate->z_) {}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
double Coordinate::GetX()
{
    return x_;
}

double Coordinate::GetY()
{
    return y_;
}

double Coordinate::GetZ()
{
    return z_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
void Coordinate::SetX(double x)
{
    x_ = x;
}

void Coordinate::SetY(double y)
{
    y_ = y;
}

void Coordinate::SetZ(double z)
{
    z_ = z;
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void Coordinate::Translate(double x, double y, double z)
{
    x_ += x;
    y_ += y;
    z_ += z;
}

bool Coordinate::CompareTo(Coordinate coordinate)
{
    if(x_ == coordinate.x_ && y_ == coordinate.y_ && z_ == coordinate.z_)
        return true;
    else
        return false;
}

double Coordinate::Distance(Coordinate coordinate)
{
    double dist = (x_ - coordinate.x_) * (x_ - coordinate.x_) + (y_ - coordinate.y_) * (y_ - coordinate.y_) + (z_ - coordinate.z_) * (z_ - coordinate.z_);
    return sqrt(dist);
}

double Coordinate::length()
{
    double length = (x_ * x_) + (y_ * y_) + (z_ * z_);
    return sqrt(length);
}
void Coordinate::Normalize()
{
    double length = this->length();
    if(length != 0.0)
    {
        x_ = x_ / length;
        y_ = y_ / length;
        z_ = z_ / length;
    }
}
double Coordinate::DotProduct(Coordinate coordinate)
{
    return ((x_ * coordinate.x_) + (y_ * coordinate.y_) + (z_ * coordinate.z_));
}
void Coordinate::CrossProduct(Coordinate coordinate)
{
    double x = x_;
    double y = y_;
    double z = z_;
    x_ = (y * coordinate.z_) - (coordinate.y_ * z);
    y_ = (z * coordinate.x_) - (coordinate.z_ * x);
    z_ = (x * coordinate.y_) - (coordinate.x_ * y);
}
void Coordinate::operator+(Coordinate coordinate)
{
    x_ += coordinate.x_;
    y_ += coordinate.y_;
    z_ += coordinate.z_;
}
void Coordinate::operator +(double addition)
{
    x_ += addition;
    y_ += addition;
    z_ += addition;
}
void Coordinate::operator-(Coordinate coordinate)
{
    x_ -= coordinate.x_;
    y_ -= coordinate.y_;
    z_ -= coordinate.z_;
}
void Coordinate::operator /(Coordinate coordinate)
{
    x_ /= coordinate.x_;
    y_ /= coordinate.y_;
    z_ /= coordinate.z_;
}
void Coordinate::operator/(double divisor)
{
    x_ /= divisor;
    y_ /= divisor;
    z_ /= divisor;
}
void Coordinate::operator *(double multiplier)
{
    x_ *= multiplier;
    y_ *= multiplier;
    z_ *= multiplier;
}
void Coordinate::TranslateAll(CoordinateVector coordinate_set, double margin, int pos)
{
    if(coordinate_set.size() == 0)
        return;
    Coordinate* direction;
    if(pos == 1)
        direction  = new Coordinate(coordinate_set.at(0));
    if(pos == -1)
        direction  = new Coordinate(coordinate_set.at(coordinate_set.size()-1));
    direction->operator -(this);
    Coordinate* offset = new Coordinate(direction);
    offset->Normalize();
    offset->operator *(margin);
    offset->operator -(direction);
    for(CoordinateVector::iterator it = coordinate_set.begin(); it != coordinate_set.end(); it++)
        (*it)->Translate(offset->GetX(), offset->GetY(), offset->GetZ());
}

void Coordinate::RotateAngularAll(CoordinateVector coordinate_set, double angle, int pos){
    if(coordinate_set.size() < 2)
        return;
    double current_angle = 0.0;
    Coordinate* a1 = this;
    Coordinate* a2;
    Coordinate* a3;
    if(pos == 1){
        a2 = new Coordinate(coordinate_set.at(0));
        a3 = new Coordinate(coordinate_set.at(1));
    }
    if(pos == -1){
        a2 = new Coordinate(coordinate_set.at(coordinate_set.size()-1));
        a3 = new Coordinate(coordinate_set.at(coordinate_set.size()-2));
    }

    Coordinate* b1 = new Coordinate(*a1);
    b1->operator -(*a2);
    Coordinate* b2 = new Coordinate(*a3);
    b2->operator -(*a2);

    current_angle = acos((b1->DotProduct(*b2)) / (b1->length() * b2->length() + DIST_EPSILON));
    double rotation_angle = ConvertDegree2Radian(angle) - current_angle;

    Coordinate* direction = new Coordinate(*b1);
    direction->CrossProduct(*b2);
    direction->Normalize();
    double** rotation_matrix = GenerateRotationMatrix(direction, a2, rotation_angle);

    if(pos == 1)
    {
        for(CoordinateVector::iterator it = coordinate_set.begin() + 1; it != coordinate_set.end(); it++)
        {
            Coordinate* coordinate = (*it);
            Coordinate* result = new Coordinate();
            result->SetX(rotation_matrix[0][0] * coordinate->GetX() + rotation_matrix[0][1] * coordinate->GetY() +
                    rotation_matrix[0][2] * coordinate->GetZ() + rotation_matrix[0][3]);
            result->SetY(rotation_matrix[1][0] * coordinate->GetX() + rotation_matrix[1][1] * coordinate->GetY() +
                    rotation_matrix[1][2] * coordinate->GetZ() + rotation_matrix[1][3]);
            result->SetZ(rotation_matrix[2][0] * coordinate->GetX() + rotation_matrix[2][1] * coordinate->GetY() +
                    rotation_matrix[2][2] * coordinate->GetZ() + rotation_matrix[2][3]);

            (*it)->SetX(result->GetX());
            (*it)->SetY(result->GetY());
            (*it)->SetZ(result->GetZ());
        }
    }
    if(pos == -1)
    {
        for(CoordinateVector::iterator it = coordinate_set.begin(); it != coordinate_set.end() - 1; it++)
        {
            Coordinate* coordinate = (*it);
            Coordinate* result = new Coordinate();
            result->SetX(rotation_matrix[0][0] * coordinate->GetX() + rotation_matrix[0][1] * coordinate->GetY() +
                    rotation_matrix[0][2] * coordinate->GetZ() + rotation_matrix[0][3]);
            result->SetY(rotation_matrix[1][0] * coordinate->GetX() + rotation_matrix[1][1] * coordinate->GetY() +
                    rotation_matrix[1][2] * coordinate->GetZ() + rotation_matrix[1][3]);
            result->SetZ(rotation_matrix[2][0] * coordinate->GetX() + rotation_matrix[2][1] * coordinate->GetY() +
                    rotation_matrix[2][2] * coordinate->GetZ() + rotation_matrix[2][3]);

            (*it)->SetX(result->GetX());
            (*it)->SetY(result->GetY());
            (*it)->SetZ(result->GetZ());
        }
    }
}

void Coordinate::RotateTorsionalAll(CoordinateVector coordinate_set, double torsion, int pos)
{
    if(coordinate_set.size() < 3)
        return;
    double current_dihedral = 0.0;
    Coordinate* a1 = this;
    Coordinate* a2;
    Coordinate* a3;
    Coordinate* a4;
    if(pos == 1)
    {
        a2 = new Coordinate(coordinate_set.at(0));
        a3 = new Coordinate(coordinate_set.at(1));
        a4 = new Coordinate(coordinate_set.at(2));
    }
    if(pos == -1)
    {
        a2 = new Coordinate(coordinate_set.at(coordinate_set.size()-1));
        a3 = new Coordinate(coordinate_set.at(coordinate_set.size()-2));
        a4 = new Coordinate(coordinate_set.at(coordinate_set.size()-3));
    }

    Coordinate* b1 = new Coordinate(*a2);
    b1->operator -(*a1);
    Coordinate* b2 = new Coordinate(*a3);
    b2->operator -(*a2);
    Coordinate* b3 = new Coordinate(*a4);
    b3->operator -(*a3);
    Coordinate* b4 = new Coordinate(*b2);
    b4->operator *(-1);

    Coordinate* b2xb3 = new Coordinate(*b2);
    b2xb3->CrossProduct(*b3);

    Coordinate* b1_m_b2n = new Coordinate(*b1);
    b1_m_b2n->operator *(b2->length());

    Coordinate* b1xb2 = new Coordinate(*b1);
    b1xb2->CrossProduct(*b2);

    current_dihedral = atan2(b1_m_b2n->DotProduct(*b2xb3), b1xb2->DotProduct(*b2xb3));

    double** torsion_matrix = GenerateRotationMatrix(b4, a2, current_dihedral - ConvertDegree2Radian(torsion));


    if(pos == 1)
    {
        for(CoordinateVector::iterator it = coordinate_set.begin() + 2; it != coordinate_set.end(); it++)
        {
            Coordinate* coordinate = (*it);
            Coordinate* result = new Coordinate();
            result->SetX(torsion_matrix[0][0] * coordinate->GetX() + torsion_matrix[0][1] * coordinate->GetY() +
                    torsion_matrix[0][2] * coordinate->GetZ() + torsion_matrix[0][3]);
            result->SetY(torsion_matrix[1][0] * coordinate->GetX() + torsion_matrix[1][1] * coordinate->GetY() +
                    torsion_matrix[1][2] * coordinate->GetZ() + torsion_matrix[1][3]);
            result->SetZ(torsion_matrix[2][0] * coordinate->GetX() + torsion_matrix[2][1] * coordinate->GetY() +
                    torsion_matrix[2][2] * coordinate->GetZ() + torsion_matrix[2][3]);

            (*it)->SetX(result->GetX());
            (*it)->SetY(result->GetY());
            (*it)->SetZ(result->GetZ());
        }
    }
    if(pos == -1)
    {
        for(CoordinateVector::iterator it = coordinate_set.begin(); it != coordinate_set.end() - 2; it++)
        {
            Coordinate* coordinate = (*it);
            Coordinate* result = new Coordinate();
            result->SetX(torsion_matrix[0][0] * coordinate->GetX() + torsion_matrix[0][1] * coordinate->GetY() +
                    torsion_matrix[0][2] * coordinate->GetZ() + torsion_matrix[0][3]);
            result->SetY(torsion_matrix[1][0] * coordinate->GetX() + torsion_matrix[1][1] * coordinate->GetY() +
                    torsion_matrix[1][2] * coordinate->GetZ() + torsion_matrix[1][3]);
            result->SetZ(torsion_matrix[2][0] * coordinate->GetX() + torsion_matrix[2][1] * coordinate->GetY() +
                    torsion_matrix[2][2] * coordinate->GetZ() + torsion_matrix[2][3]);

            (*it)->SetX(result->GetX());
            (*it)->SetY(result->GetY());
            (*it)->SetZ(result->GetZ());
        }
    }
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void Coordinate::Print(std::ostream& out)
{
    if(this->CompareTo(Coordinate(dNotSet, dNotSet, dNotSet)) == true)
        out << std::setw(10) << " " << ", " << std::setw(10) << " " << ", " << std::setw(10) << " ";
    else
        out << std::setw(10) << x_ << ", " << std::setw(10) << y_ << ", " << std::setw(10) << z_;
}


