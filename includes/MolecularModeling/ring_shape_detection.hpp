#ifndef ring_shape_detection_HPP
#define ring_shape_detection_HPP

#include <stdio.h>
#include <sstream>
#include <string.h>
#include "../GeometryTopology/coordinate.hpp"
#include "../Glycan/monosaccharide.hpp"
#include "../common.hpp"

namespace glylib
{
    //****************************************************//
    //*************** Glylib structs *********************//
    typedef struct
    {
        char *N; /// the name of the file -- use strdup to allocate/copy
        FILE *F; /// the file pointer
    } fileset;
    typedef struct
    {
        double A,B,C,D;
    } plane; ///< Standard cartesian plane with equation Ax+By+Cz+D=0
    typedef struct
    {
        double i,j,k,d; ///< i,j,k vector with magnitude d
    } vectormag_3D; ///< Vector manipulations often require magnitudes, so this struct contains that information.
    //****************************************************//
}

//****************************************************//
//*************** Glylib functions ******************//
namespace glylib
{
    typedef std::vector<plane> PlaneVector;
    double get_signed_distance_from_point_to_plane(plane p, GeometryTopology::Coordinate pt);
    vectormag_3D get_crossprod(vectormag_3D a, vectormag_3D b);
    std::string calculateRingShapeBFMP(GeometryTopology::CoordinateVector ring_coordinates, int cut_off = 10);
    double calculateTorsionAngle(GeometryTopology::Coordinate *coord1, GeometryTopology::Coordinate *coord2, GeometryTopology::Coordinate *coord3, GeometryTopology::Coordinate *coord4);
    plane get_plane_for_ring(int n, GeometryTopology::CoordinateVector r);
    FILE *myfopen(const char* Name,const char* option); ///< Opens a file with error checking.
    //****************************************************//

    //****************************************************//
    //*** Oliver functions to reduce repeition (a bit) ***//
    bool does_conformation_match(int jloops, int *sortedplanes, int *atoms, int atoms_size, int *list, int required_matches);
    int count_conformation_matches(int jloops, int *sortedplanes, int *atoms, int atoms_size, int *list);
    //****************************************************//

    //****************************************************//
    //************* The BFMP function ********************//

    std::string CalculateRingShapeBFMP( Glycan::Monosaccharide* mono );
    std::string CalculateRingShapeBFMP(GeometryTopology::CoordinateVector ring_coordinates, int cut_off = 10);
    //****************************************************//
}

// This doesn't belong here. Doesn't belong in Assembly (where I took it from), it belongs in Monosaccharide right? Except that's not a class.
namespace gmml
{
    GeometryTopology::CoordinateVector GetCycleAtomCoordinates( Glycan::Monosaccharide* mono );
}
#endif // ring_shape_detection_HPP
