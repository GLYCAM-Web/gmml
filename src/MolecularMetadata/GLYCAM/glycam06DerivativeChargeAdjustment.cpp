#include "../../../includes/MolecularMetadata/GLYCAM/glycam06DerivativeChargeAdjustment.hpp"
#include <iostream> // for cout, can remove after debug
#include <locale> // for isLower()
#include <sstream> // for stringstream

using gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeChargeAdjustmentLookupContainer;

//////////////////////////////////////////////////////////
//                      QUERY FUNCTIONS                 //
//////////////////////////////////////////////////////////
    std::string glycamResidueCode_ ;
    std::string adjustmentAtom_ ;
    double charge_ ;
double Glycam06DerivativeChargeAdjustmentLookupContainer::GetAdjustmentCharge(std::string queryResidueName)
{ // e.g. input = SO3, output = +0.008.
    for (const auto& entry : glycam06DerivativeChargeAdjustmentLookup_)
    {
        if (entry.glycamResidueCode_ == queryResidueName)
            return entry.charge_;
    }
    return 0.0;
}

std::string Glycam06DerivativeChargeAdjustmentLookupContainer::GetAdjustmentAtom(std::string queryResidueName)
{   // e.g. input = SO3, output = O.
    for (const auto& entry : glycam06DerivativeChargeAdjustmentLookup_)
    {
        if (entry.glycamResidueCode_ == queryResidueName)
            return entry.adjustmentAtom_;
    }
    return "O";
}

//////////////////////////////////////////////////////////
//                    INITIALIZER                       //
//////////////////////////////////////////////////////////

Glycam06DerivativeChargeAdjustmentLookupContainer::Glycam06DerivativeChargeAdjustmentLookupContainer()
{
    glycam06DerivativeChargeAdjustmentLookup_ =
    {// residueCode_, adjustmentAtom_ , charge_
        {"ACX", "C",  0.008},
        {"MEX", "C", -0.039},
        {"SO3", "O", +0.031},
    };
}
