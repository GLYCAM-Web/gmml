
#include "../../../includes/FileSet/TopologyFileSpace/topologyatomtype.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAtomType::TopologyAtomType() {}

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
    coefficient_a_ = coefficient_a;
}
void TopologyAtomType::SetCoefficientB(TopologyCoefficientBMap coefficient_b)
{
    coefficient_b_ = coefficient_b;
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

