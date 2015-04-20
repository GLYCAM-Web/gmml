#ifndef COMMON_HPP
#define COMMON_HPP

#include <string>
#include <vector>
#include <math.h>
#include <map>

#include "Geometry/coordinate.hpp"

namespace gmml
{

    //*******************************************
    typedef std::map<std::string, std::string> ResidueNameMap;
    typedef std::map<std::string, std::vector<std::string> > ResidueNameAtomNamesMap;
    typedef Geometry::Coordinate Vector;

    //*******************************************

    const double dNotSet = 123456789.0;
    const int iNotSet = -123456;
    const int iPdbLineLength = 80;
    const double dSulfurCutoff = 2.5;
    const double PI_RADIAN = 4.0*atan(1.0);
    const double PI_DEGREE = 180.0;
    const double EPSILON = 0.001;
    const double dCutOff = 1.6;
    const int PdbResidueThreshold = 500;
    const int DEFAULT_DUMMY_ATOMS = 3;
    const char BLANK_SPACE = '?';

    /*! \enum
      * Enumerator to possible n chain termination
      */
    enum PossibleNChainTermination
    {
        COCH3 = 1,
        NH3 = 2
    };

    /*! \enum
      * Enumerator to possible c chain termination
      */
    enum PossibleCChainTermination
    {
        NH2 = 1,
        NHCH3 = 2,
        CO2 = 3
    };

    /*! \enum
      * Enumerator to possible HIS mapping in pdb preprocessor
      */
    enum PdbPreprocessorHISMapping
    {
        HIE = 1,
        HIP = 2,
        HID = 3
    };

    /*! \enum
      * Enumerator to possible input file type to central data structure
      */
    enum InputFileType
    {
        PDB,
        LIB,
        PREP,
        TOP,
        TOP_CRD,
        MULTIPLE,
        UNKNOWN
    };

    /*! \enum
      * Enumerator to possible options for building a graph for a central data structure
      */
    enum BuildingStructureOption
    {
        DISTANCE,
        ORIGINAL,
        DATABASE
    };

    /*! \enum
      * Enumerator to parameter file type
      */
    enum ParameterFileType
    {
        MAIN = 0,           /*!< Main parameter file >*/
        MODIFIED = 1        /*!< Force modified parameter file >*/
    };

    /*! \enum
      * Topological type enumerator
      */
    enum TopologicalType
    {
        kTopTypeE,
        kTopTypeS,
        kTopTypeB,
        kTopType3,
        kTopType4,
        kTopTypeM
    };

    /*! \enum
      * Enumerator to topological type stack
      */
    enum TopologicalTypeStackElement
    {
        EMPTY = 0,
        M1 = 1,
        M2 = 2,
        M3 = 3,
        M4 = 4,
        S1 = 5,
        B1 = 6,
        B2 = 7,
        T1 = 8,
        T2 = 9,
        T3 = 10
    };

    /*! \enum
      * Enumerator to topological type stack
      */
    enum GraphSearchNodeStatus
    {
        UNVISITED = 0,
        VISITED = 1,
        DONE = 2
    };

}

#endif // COMMON_HPP
