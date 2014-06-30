
#include "../../../includes/FileSet/TopologyFileSpace/topologyatompair.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAtomPair::TopologyAtomPair() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int TopologyAtomPair::GetAtomTypeIndex()
{
    return atom_type_index_;
}
string TopologyAtomPair::GetAtomType()
{
    return atom_type_;
}
TopologyAtomPair::TopologyCoefficientMap TopologyAtomPair::GetCoefficientA()
{
    return coefficient_a_;
}
TopologyAtomPair::TopologyCoefficientMap TopologyAtomPair::GetCoefficientB()
{
    return coefficient_b_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAtomPair::SetAtomTypeIndex(int atom_type_index)
{
    atom_type_index_ = atom_type_index;
}
void TopologyAtomPair::SetAtomType(string atom_type)
{
    atom_type_ = atom_type;
}
void TopologyAtomPair::SetCoefficientA(TopologyCoefficientMap coefficient_a)
{
    coefficient_a_.clear();
    for(TopologyCoefficientMap::iterator it = coefficient_a.begin(); it != coefficient_a.end(); it++)
    {
       double coefficienta = (*it).second;
       int atom_type_index = (*it).first;
       coefficient_a_[atom_type_index] = coefficienta;
    }
}
void TopologyAtomPair::SetCoefficientB(TopologyCoefficientMap coefficient_b)
{
    coefficient_b_.clear();
    for(TopologyCoefficientMap::iterator it = coefficient_b.begin(); it != coefficient_b.end(); it++)
    {
       double coefficientb = (*it).second;
       int atom_type_index = (*it).first;
       coefficient_b_[atom_type_index] = coefficientb;
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAtomPair::Print(ostream &out)
{
}

