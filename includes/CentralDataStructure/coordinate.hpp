#ifndef INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP

#include <iostream>

namespace cds
{
class Coordinate
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    Coordinate();
    Coordinate(double x, double y, double z);
    Coordinate(const Coordinate& c);
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const double& getX() const {return x_;}
    inline const double& getY() const {return y_;}
    inline const double& getZ() const {return z_;}
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setX(double x) {x_ = x;}
    inline void setY(double y) {y_ = y;}
    inline void setZ(double z) {z_ = z;}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    void translate(const double x, const double y, const double z);
    double distance(const Coordinate &c) const;
    double getLength() const;
    void normalize();
    double dotProduct(const Coordinate& c);
    void crossProduct(const Coordinate& c);
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
    void Print(std::ostream& out = std::cerr);
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    double x_;          /*!< x */
    double y_;          /*!< y */
    double z_;          /*!< z */
};
}

#endif // INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP
