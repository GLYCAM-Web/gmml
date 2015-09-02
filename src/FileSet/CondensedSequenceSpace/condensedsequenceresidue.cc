#include "../../../includes/FileSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace CondensedSequenceResidueSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequenceResidue::CondensedSequenceResidue()
{
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
bool CondensedSequenceResidue::GetIsTerminal()
{
    return is_terminal_;
}
string CondensedSequenceResidue::GetIsomer()
{
    return isomer_;
}
string CondensedSequenceResidue::GetConfiguration()
{
    return configuration_;
}
string CondensedSequenceResidue::GetName()
{
    return name_;
}
int CondensedSequenceResidue::GetAnomericCarbon()
{
    return anomeric_carbon_;
}
int CondensedSequenceResidue::GetOxygenPosition()
{
    return oxygen_position_;
}
CondensedSequenceResidue::DerivativeMap CondensedSequenceResidue::GetDerivatives()
{
    return derivatives_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequenceResidue::SetIsTerminal(bool is_terminal)
{
    is_terminal_ = is_terminal;
}
void CondensedSequenceResidue::SetIsomer(string isomer)
{
    isomer_ = isomer;
}
void CondensedSequenceResidue::SetConfiguration(string configuration)
{
    configuration_ = configuration;
}
void CondensedSequenceResidue::SetName(string name)
{
    name_ = name;
}
void CondensedSequenceResidue::SetAnomericCarbon(int anomeric_carbon)
{
    anomeric_carbon_ = anomeric_carbon;
}
void CondensedSequenceResidue::SetDerivatives(DerivativeMap derivatives)
{
    derivatives_.clear();
    for(DerivativeMap::iterator it = derivatives.begin(); it != derivatives.end(); it++)
    {
        int index = (*it).first;
        string derivative = (*it).second;
        derivatives_[index] = derivative;
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequenceResidue::Print(ostream &out)
{}





