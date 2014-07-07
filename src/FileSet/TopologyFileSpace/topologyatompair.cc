
#include "../../../includes/FileSet/TopologyFileSpace/topologyatompair.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAtomPair::TopologyAtomPair() {}

TopologyAtomPair::TopologyAtomPair(string pair_type, double coefficient_a, double coefficient_b) :
    pair_type_(pair_type), coefficient_a_(coefficient_a), coefficient_b_(coefficient_b) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string TopologyAtomPair::GetPairType()
{
    return pair_type_;
}
double TopologyAtomPair::GetCoefficientA()
{
    return coefficient_a_;
}
double TopologyAtomPair::GetCoefficientB()
{
    return coefficient_b_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAtomPair::SetPairType(string pair_type)
{
    pair_type_ = pair_type;
}
void TopologyAtomPair::SetCoefficientA(double coefficient_a)
{
    coefficient_a_ = coefficient_a;
}
void TopologyAtomPair::SetCoefficientB(double coefficient_b)
{
    coefficient_b_ = coefficient_b;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAtomPair::Print(ostream &out)
{
    out << pair_type_ << ": a = " << coefficient_a_ << "; b = " << coefficient_b_ << endl;
}

