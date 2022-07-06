#ifndef GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06DERIVATIVECHARGEADJUSTMENT_HPP
#define GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06DERIVATIVECHARGEADJUSTMENT_HPP
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
    std::string glycamResidueCode_ ;
    std::string adjustmentAtom_ ;
    double charge_ ;
} ;

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
} // close namespace
} // close namespace
} // close namespace
#endif // GMML_INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06DERIVATIVECHARGEADJUSTMENT_HPP
