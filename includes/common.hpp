#ifndef COMMON_HPP
#define COMMON_HPP

#include <string>
#include <vector>
#include <math.h>

namespace gmml
{
    const double dNotSet = 123456789.0;
    const int iNotSet = -123456;
    const int iPdbLineLength = 80;
    const double dSulfurCutoff = 2.5;
    const double PI_RADIAN = 4.0*atan(1.0);
    const double PI_DEGREE = 180.0;
    const double EPSILON = 0.001;
    const double dCutOff = 1.6;

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

}

#endif // COMMON_HPP
