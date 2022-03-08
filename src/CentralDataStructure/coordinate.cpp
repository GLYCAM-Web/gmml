#include <cmath> //sqrt

#include "includes/CentralDataStructure/coordinate.hpp"

using cds::Coordinate;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
Coordinate::Coordinate() : x_(0.0), y_(0.0), z_(0.0) {}

Coordinate::Coordinate(double x, double y, double z) : x_(x), y_(y), z_(z) {}

Coordinate::Coordinate(const Coordinate &c) : x_(c.getX()), y_(c.getY()), z_(c.getZ()) {}
//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////
void Coordinate::translate(const double x, const double y, const double z)
{
    x_ += x;
    y_ += y;
    z_ += z;
    return;
}

double Coordinate::distance(const Coordinate &c)const
{
    double dist = (x_ - c.getX()) * (x_ - c.getX()) + (y_ - c.getY()) * (y_ - c.getY()) + (z_ - c.getZ()) * (z_ - c.getZ());
    return sqrt(dist);
}

double Coordinate::getLength() const
{
    return sqrt((x_ * x_) + (y_ * y_) + (z_ * z_));
}

void Coordinate::normalize()
{
    double length = this->getLength();
    if(length != 0.0)
    {
        x_ = x_ / length;
        y_ = y_ / length;
        z_ = z_ / length;
    }
    return;
}

double Coordinate::dotProduct(const Coordinate& c)
{
    return ((x_ * c.getX()) + (y_ * c.getY()) + (z_ * c.getZ()));
}

void Coordinate::crossProduct(const Coordinate& c)
{
    x_ = (y_ * c.getZ()) - (c.getY() * z_);
    y_ = (z_ * c.getX()) - (c.getZ() * x_);
    z_ = (x_ * c.getY()) - (c.getX() * y_);
    return;
}
//////////////////////////////////////////////////////////
//                    DISPLAY                           //
//////////////////////////////////////////////////////////
void Coordinate::Print(std::ostream& out)
{
    out << "" << this->getX() << "  " << this->getY() << "  " << this->getX() << "\n";
    return;
}
