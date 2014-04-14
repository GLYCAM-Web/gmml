#include <math.h>

#include "../includes/Geometry/coordinate.hpp"
#include "../includes/common.hpp"

using namespace Geometry;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
Coordinate::Coordinate() : x_(0.0), y_(0.0), z_(0.0) {}

Coordinate::Coordinate(double x, double y, double z) : x_(x), y_(y), z_(z) {}

Coordinate::Coordinate(const Coordinate &coordinate) : x_(coordinate.x_), y_(coordinate.y_), z_(coordinate.z_) {}

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


