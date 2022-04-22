//#ifndef INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP
//#define INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP
//
//#include <iostream>
//
//namespace cds
//{
//class cdsCoordinate
//{
//public:
//    //////////////////////////////////////////////////////////
//    //                    CONSTRUCTOR                       //
//    //////////////////////////////////////////////////////////
//    cdsCoordinate();
//    cdsCoordinate(double x, double y, double z);
//    cdsCoordinate(const cdsCoordinate& c);
////    //////////////////////////////////////////////////////////
////    //                    OPERATOR                          //
////    //////////////////////////////////////////////////////////
////    cdsCoordinate& operator= (const cdsCoordinate& c);
////    cdsCoordinate* operator= (const cdsCoordinate& c);
//    //////////////////////////////////////////////////////////
//    //                    ACCESSOR                          //
//    //////////////////////////////////////////////////////////
//    inline const double& getX() const {return x_;}
//    inline const double& getY() const {return y_;}
//    inline const double& getZ() const {return z_;}
//    //////////////////////////////////////////////////////////
//    //                    MUTATOR                           //
//    //////////////////////////////////////////////////////////
//    inline void setX(double x) {x_ = x;}
//    inline void setY(double y) {y_ = y;}
//    inline void setZ(double z) {z_ = z;}
//    //////////////////////////////////////////////////////////
//    //                    FUNCTIONS                         //
//    //////////////////////////////////////////////////////////
//    void translate(const double x, const double y, const double z);
//    double distance(const cdsCoordinate &c) const;
//    double getLength() const;
//    void normalize();
//    double dotProduct(const cdsCoordinate& c);
//    void crossProduct(const cdsCoordinate& c);
//    //////////////////////////////////////////////////////////
//    //                    DISPLAY                           //
//    //////////////////////////////////////////////////////////
//    void Print(std::ostream& out = std::cerr);
//private:
//    //////////////////////////////////////////////////////////
//    //                    ATTRIBUTES                        //
//    //////////////////////////////////////////////////////////
//    double x_;          /*!< x */
//    double y_;          /*!< y */
//    double z_;          /*!< z */
//};
//}
//
//#endif // INCLUDES_CENTRALDATASTRUCTURE_COORDINATE_HPP
