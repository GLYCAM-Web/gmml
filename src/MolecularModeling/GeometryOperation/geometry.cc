#include "../../../includes/MolecularModeling/GeometryOperation/geometry.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/GeometryTopology/coordinate.hpp"

#include <iostream>
#include <math.h>

using namespace GeometryOperation;
using namespace GeometryTopology;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Geometry::Geometry(){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
Geometry::CoordinateVector Geometry::RotateCoordinates(Coordinate* pivot_point,Coordinate* direction_point,double rotation_angle,CoordinateVector coordinate_set)
{
     //Please refer https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions for the technique used in coordinate rotation with respect to arbitrary axis.

    //direction point/vector point of Arbitrary Axis
    double u = direction_point->GetX();
    double v = direction_point->GetY();
    double w = direction_point->GetZ();

    //Check whether a vector's length is less than 0
    double l=sqrt(u*u+v*v+w*w);
    if(l<0)
    {
        cout<<"Direction vector length is less than 1"<<endl;
    }else
    {
        // Normalize the direction vector.
        u = u/l;
        v = v/l;
        w = w/l;


        // Set some intermediate values.
        double u2 = u*u;
        double v2 = v*v;
        double w2 = w*w;


        //start point of Arbitrary Axis
        double a = pivot_point->GetX();
        double b = pivot_point->GetY();
        double c = pivot_point->GetZ();

        //converting degree to radians
        double angle = ConvertDegree2Radian(rotation_angle);

        double cosT = cos(angle);
        double oneMinusCosT = 1 - cosT;
        double sinT = sin(angle);

        for(CoordinateVector::iterator it = coordinate_set.begin(); it != coordinate_set.end(); it++)
        {
            Coordinate* coordinate = (*it);

            double x=coordinate->GetX();
            double y=coordinate->GetY();
            double z= coordinate->GetZ();

            Coordinate* result = new Coordinate();
            result->SetX((a*(v2 + w2) - u*(b*v + c*w - u*x - v*y - w*z)) * oneMinusCosT+ x*cosT+ (-c*v + b*w - w*y + v*z)*sinT);
            result->SetY((b*(u2 + w2) - v*(a*u + c*w - u*x - v*y - w*z)) * oneMinusCosT+ y*cosT+ (c*u - a*w + w*x - u*z)*sinT);
            result->SetZ((c*(u2 + v2) - w*(a*u + b*v - u*x - v*y - w*z)) * oneMinusCosT+ z*cosT+ (-b*u + a*v - v*x + u*y)*sinT);

            (*it)->SetX(result->GetX());
            (*it)->SetY(result->GetY());
            (*it)->SetZ(result->GetZ());
        }

           return coordinate_set;
    }

}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
