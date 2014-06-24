#ifndef COMMON_HPP
#define COMMON_HPP

#include <string>
#include <vector>

namespace gmml
{
    const double dNotSet = 123456789.0;
    const int iNotSet = -123456;
    const int iPdbLineLength = 80;
    const double dSulfurCutoff = 2.5;

    enum PossibleNChainTermination
    {
        COCH3 = 1,
        NH3 = 2
    };

    enum PossibleCChainTermination
    {
        NH2 = 1,
        NHCH3 = 2,
        CO2 = 3
    };

    enum PdbPreprocessorHISMapping
    {
        HIE = 1,
        HIP = 2,
        HID = 3
    };

    enum InputFileType
    {
        PDB,
        LIB,
        PREP,
        PARAM,
        TOP
    };

    enum GraphType
    {
        RESIDUE,
        MOLECULE,
        ASSEMBLY
    };

}

#endif // COMMON_HPP
