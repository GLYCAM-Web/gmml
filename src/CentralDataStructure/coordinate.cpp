#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <cmath>
#include <limits>
#include <iomanip> // setw

using cds::Coordinate;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
Coordinate::Coordinate(double x, double y, double z) : x_(x), y_(y), z_(z) {}
Coordinate::Coordinate(const std::string x, const std::string y, const std::string z)
{
    try
    {
        x_ = std::stod(x);
        y_ = std::stod(y);
        z_ = std::stod(z);
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Could not convert these strings to doubles: " + x + ", " + y + ", " + z + ", ");
        throw;
    }
}
Coordinate::Coordinate(const Coordinate &coordinate) : x_(coordinate.x_), y_(coordinate.y_), z_(coordinate.z_) {}
Coordinate::Coordinate(Coordinate* coordinate) : x_(coordinate->x_), y_(coordinate->y_), z_(coordinate->z_) {}
//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void Coordinate::Translate(double x, double y, double z)
{
    x_ += x;
    y_ += y;
    z_ += z;
}

bool Coordinate::withinDistance(const Coordinate *coordinate, const double distance) const
{ // Vast majority of calls to this function will be able to return false after first if.
    if (this->GetX() - coordinate->GetX() < distance)
    {
        if (this->GetY() - coordinate->GetY() < distance)
        {
            if (this->Distance(*coordinate) < distance) // Need to test if also checking Z is faster.
            {
                return true;
            }
        }
    }
    return false;
}

double Coordinate::Distance(const Coordinate &coordinate) const
{
    double dist = (x_ - coordinate.x_) * (x_ - coordinate.x_) + (y_ - coordinate.y_) * (y_ - coordinate.y_) + (z_ - coordinate.z_) * (z_ - coordinate.z_);
    if(dist > 0.00000001) // can sometimes measure distance to self, in which case get sqrt(0), which should be fine but zero is funky and somtimes is actually slightly negative.
    {
        return sqrt(dist);
    }
    return 0.0;
}

double Coordinate::length() const
{
    return sqrt( (this->GetX() * this->GetX()) + (this->GetY() * this->GetY()) + (this->GetZ() * this->GetZ()) );
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
// Should this not return a coord instead of altering this one?
void Coordinate::CrossProduct(Coordinate coordinate)
{
    double x = x_;
    double y = y_;
    double z = z_;
    x_ = (y * coordinate.z_) - (coordinate.y_ * z);
    y_ = (z * coordinate.x_) - (coordinate.z_ * x);
    z_ = (x * coordinate.y_) - (coordinate.x_ * y);
}
// These can all be better, the call is weird:
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
//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void Coordinate::Print(std::ostream& out) const
{
    if(this->GetX() == constants::dNotSet || this->GetY() == constants::dNotSet || this->GetZ() == constants::dNotSet )
        out << std::setw(10) << " " << ", " << std::setw(10) << " " << ", " << std::setw(10) << " ";
    else
        out << std::setw(10) << x_ << ", " << std::setw(10) << y_ << ", " << std::setw(10) << z_;
}

std::string Coordinate::ToString() const
{
    std::stringstream ss;
    this->Print(ss);
    return ss.str();
}

