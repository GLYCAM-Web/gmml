#ifndef GLYCAM06_DERIVATIVE_CHARGE_ADJUSTMENT_HPP
#define GLYCAM06_DERIVATIVE_CHARGE_ADJUSTMENT_HPP
#include <string>
#include <vector>

namespace gmml
{
    namespace MolecularMetadata
    {
        namespace GLYCAM
        {
            struct ResidueChargeAdjustmentData
            {
                std::string glycamResidueCode_;
                std::string adjustmentAtom_;
                double charge_;
            };

            class Glycam06DerivativeChargeAdjustmentLookupContainer
            {
              public:
                //////////////////////////////////////////////////////////
                //                       CONSTRUCTOR                    //
                //////////////////////////////////////////////////////////
                Glycam06DerivativeChargeAdjustmentLookupContainer();
                //////////////////////////////////////////////////////////
                //                      QUERY FUNCTIONS                 //
                //////////////////////////////////////////////////////////
                double GetAdjustmentCharge(std::string queryResidueCode);
                std::string GetAdjustmentAtom(std::string queryResidueCode);

              private:
                std::vector<ResidueChargeAdjustmentData> glycam06DerivativeChargeAdjustmentLookup_;
            };
        } // namespace GLYCAM
    }     // namespace MolecularMetadata
} // namespace gmml
#endif // GLYCAM06_DERIVATIVE_CHARGE_ADJUSTMENT_HPP
