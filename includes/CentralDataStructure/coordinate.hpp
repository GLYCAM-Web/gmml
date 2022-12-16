#ifndef INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP

#include "includes/CodeUtils/constants.hpp"

#include <iostream>
#include <vector>

namespace cds
{
class Coordinate
{
public:
    //////////////////////////////////////////////////////////
    //                       Constructor                    //
    //////////////////////////////////////////////////////////
    Coordinate() {}
    Coordinate(const std::string x, const std::string y, const std::string z);
    Coordinate(double x, double y, double z);
    //Coordinate(const Coordinate& coordinate);
    //Coordinate(Coordinate* coordinate);
    //////////////////////////////////////////////////////////
    //                           ACCESSOR                   //
    //////////////////////////////////////////////////////////
    double GetX() const {return x_;}
    double GetY() const {return y_;}
    double GetZ() const {return z_;}
    //////////////////////////////////////////////////////////
    //                           MUTATOR                    //
    //////////////////////////////////////////////////////////
    void SetX(const double x) {x_ = x;}
    void SetY(const double y) {y_ = y;}
    void SetZ(const double z) {z_ = z;}
    void Set(const double x, const double y, const double z) {x_ = x; y_ = y; z_ = z;}
    //////////////////////////////////////////////////////////
    //                         FUNCTIONS                    //
    //////////////////////////////////////////////////////////
    void Translate(const double x, const double y, const double z);
    bool withinDistance(const Coordinate *coordinate, const double distance) const;
    double Distance(const Coordinate *coordinate) const;
    double length() const;
    void Normalize();
    double DotProduct(Coordinate coordinate);
    void CrossProduct(Coordinate coordinate);
    void operator+(Coordinate coordinate);
    void operator +(double addition);
    void operator-(Coordinate coordinate);
    void operator /(Coordinate coordinate);
    void operator/(double divisor);
    void operator*(double multiplier);
    bool operator == (const Coordinate& rhs) const { return (this->GetX() == rhs.GetX() && this->GetY() == rhs.GetY() && this->GetZ() == rhs.GetZ());}
    Coordinate& operator += (const Coordinate& rhs) { x_ += rhs.GetX();
                                                            y_ += rhs.GetY();
                                                            z_ += rhs.GetZ();
                                                            return *this;
    }

    //////////////////////////////////////////////////////////
    //                     DISPLAY FUNCTIONS                //
    //////////////////////////////////////////////////////////
    void Print(std::ostream& out = std::cerr) const;
    std::string ToString() const;
private:
    //////////////////////////////////////////////////////////
    //                         ATTRIBUTES                   //
    //////////////////////////////////////////////////////////
    double x_ = constants::dNotSet;
    double y_ = constants::dNotSet;
    double z_ = constants::dNotSet;
};
} // namespace
#endif
