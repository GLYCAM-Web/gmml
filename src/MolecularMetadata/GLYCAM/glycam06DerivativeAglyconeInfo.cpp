#include "./includes/MolecularMetadata/GLYCAM/glycam06DerivativeAglyconeInfo.hpp"


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::Glycam06DerivativeAglyconeConnectionAtomLookup;

Glycam06DerivativeAglyconeConnectionAtomLookup::Glycam06DerivativeAglyconeConnectionAtomLookup()
{
    glycam06DerivativeAglyconeConnectionAtomLookup_ =
    {
        //Aglycones
        { "ROH" , "O1"  } ,
        { "OME" , "O"   } ,
        { "TBT" , "O1"  } ,
        //Aglycones amino acid
        { "NLN" , "ND2" } ,
        { "OLT" , "OG1" } ,
        { "OLS" , "OG"  } ,
        //Derivatives
        { "SO3" , "S1"  } ,
        { "MEX" , "CH3" } ,
        { "ACX" , "C2A" } ,
    }; 
}
