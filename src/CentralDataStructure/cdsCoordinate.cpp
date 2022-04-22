//#include <cmath> //sqrt
//#include "includes/CentralDataStructure/cdsCoordinate.hpp"
//
//using cds::cdsCoordinate;
////////////////////////////////////////////////////////////
////                    CONSTRUCTOR                       //
////////////////////////////////////////////////////////////
//cdsCoordinate::cdsCoordinate() : x_(0.0), y_(0.0), z_(0.0) {}
//
//cdsCoordinate::cdsCoordinate(double x, double y, double z) : x_(x), y_(y), z_(z) {}
//
//cdsCoordinate::cdsCoordinate(const cdsCoordinate &c) : x_(c.getX()), y_(c.getY()), z_(c.getZ()) {}
////////////////////////////////////////////////////////////
////                    OPERATOR                          //
////////////////////////////////////////////////////////////
////cdsCoordinate::cdsCoordinate& operator= (const cdsCoordinate* c)
////{
////
////}
////cdsCoordinate::cdsCoordinate* operator= (const cdsCoordinate& c)
////{
////
////}
////////////////////////////////////////////////////////////
////                    FUNCTIONS                         //
////////////////////////////////////////////////////////////
//void cdsCoordinate::translate(const double x, const double y, const double z)
//{
//    x_ += x;
//    y_ += y;
//    z_ += z;
//    return;
//}
//
//double cdsCoordinate::distance(const cdsCoordinate &c)const
//{
//    double dist = (x_ - c.getX()) * (x_ - c.getX()) + (y_ - c.getY()) * (y_ - c.getY()) + (z_ - c.getZ()) * (z_ - c.getZ());
//    return sqrt(dist);
//}
//
//double cdsCoordinate::getLength() const
//{
//    return sqrt((x_ * x_) + (y_ * y_) + (z_ * z_));
//}
//
//void cdsCoordinate::normalize()
//{
//    double length = this->getLength();
//    if(length != 0.0)
//    {
//        x_ = x_ / length;
//        y_ = y_ / length;
//        z_ = z_ / length;
//    }
//    return;
//}
//
//double cdsCoordinate::dotProduct(const cdsCoordinate& c)
//{
//    return ((x_ * c.getX()) + (y_ * c.getY()) + (z_ * c.getZ()));
//}
//
//void cdsCoordinate::crossProduct(const cdsCoordinate& c)
//{
//    x_ = (y_ * c.getZ()) - (c.getY() * z_);
//    y_ = (z_ * c.getX()) - (c.getZ() * x_);
//    z_ = (x_ * c.getY()) - (c.getX() * y_);
//    return;
//}
////////////////////////////////////////////////////////////
////                    DISPLAY                           //
////////////////////////////////////////////////////////////
//void cdsCoordinate::Print(std::ostream& out)
//{
//    out << "" << this->getX() << "  " << this->getY() << "  " << this->getX() << "\n";
//    return;
//}
