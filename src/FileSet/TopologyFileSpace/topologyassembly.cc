
#include "../../../includes/FileSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyresidue.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAssembly::TopologyAssembly() {}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string TopologyAssembly::GetAssemblyName()
{
    return assembly_name_;
}
TopologyAssembly::TopologyResidueMap TopologyAssembly::GetResidues()
{
    return residues_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAssembly::SetAssemblyName(const string assembly_name)
{
    assembly_name_ = assembly_name;
}
void TopologyAssembly::SetResidues(TopologyResidueMap residues)
{
    residues_.clear();
    for(TopologyResidueMap::iterator it = residues.begin(); it != residues.end(); it++)
    {
        TopologyResidue* residue = (*it).second;
        residues_[residue->GetResidueName()] = residue;
    }
}


//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAssembly::Print(ostream &out)
{
}

