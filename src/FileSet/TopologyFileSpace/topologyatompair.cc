
#include "../../../includes/FileSet/TopologyFileSpace/topologyatompair.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAtomPair::TopologyAtomPair() {}

TopologyAtomPair::TopologyAtomPair(string atom_type, TopologyCoefficientMap coefficient_a_map, TopologyCoefficientMap coefficient_b_map) :
    atom_type_(atom_type)
{
    coefficient_a_.clear();
    for(TopologyCoefficientMap::iterator it = coefficient_a_map.begin(); it != coefficient_a_map.end(); it++)
    {
        string atom_type = (*it).first;
        double coefficient = (*it).second;
        coefficient_a_[atom_type] = coefficient;
    }
    coefficient_b_.clear();
    for(TopologyCoefficientMap::iterator it = coefficient_b_map.begin(); it != coefficient_b_map.end(); it++)
    {
        string atom_type = (*it).first;
        double coefficient = (*it).second;
        coefficient_b_[atom_type] = coefficient;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
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
       string atom_type = (*it).first;
       coefficient_a_[atom_type] = coefficienta;
    }
}
void TopologyAtomPair::SetCoefficientB(TopologyCoefficientMap coefficient_b)
{
    coefficient_b_.clear();
    for(TopologyCoefficientMap::iterator it = coefficient_b.begin(); it != coefficient_b.end(); it++)
    {
       double coefficientb = (*it).second;
       string atom_type = (*it).first;
       coefficient_b_[atom_type] = coefficientb;
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

