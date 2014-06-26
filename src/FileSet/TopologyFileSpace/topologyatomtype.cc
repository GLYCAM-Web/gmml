
#include "../../../includes/FileSet/TopologyFileSpace/topologyatomtype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAtomType::TopologyAtomType() {}

TopologyAtomType::TopologyAtomType(int atom_type_index, int index) : atom_type_index_(atom_type_index), index_(index) {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
int TopologyAtomType::GetAtomTypeIndex()
{
    return atom_type_index_;
}
int TopologyAtomType::GetIndex()
{
    return index_;
}
TopologyAtomType::TopologyCoefficientAMap TopologyAtomType::GetCoefficientA()
{
    return coefficient_a_;
}
TopologyAtomType::TopologyCoefficientBMap TopologyAtomType::GetCoefficientB()
{
    return coefficient_b_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAtomType::SetAtomTypeIndex(int atom_type_index)
{
    atom_type_index_ = atom_type_index;
}
void TopologyAtomType::SetIndex(int index)
{
    index_ = index;
}
void TopologyAtomType::SetCoefficientA(TopologyCoefficientAMap coefficient_a)
{
    coefficient_a_.clear();
    for(TopologyCoefficientAMap::iterator it = coefficient_a.begin(); it != coefficient_a.end(); it++)
    {
       double coefficienta = (*it).second;
       int atom_type_index = (*it).first;
       coefficient_a_[atom_type_index] = coefficienta;
    }
}
void TopologyAtomType::SetCoefficientB(TopologyCoefficientBMap coefficient_b)
{
    coefficient_b_.clear();
    for(TopologyCoefficientBMap::iterator it = coefficient_b.begin(); it != coefficient_b.end(); it++)
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
void TopologyAtomType::Print(ostream &out)
{
}

