
#include "../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"

using namespace std;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
TopologyAssembly::TopologyAssembly()
{
    assembly_name_ = "";
    residues_ = TopologyResidueVector();
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string TopologyAssembly::GetAssemblyName()
{
    return assembly_name_;
}
TopologyAssembly::TopologyResidueVector TopologyAssembly::GetResidues()
{
    return residues_;
}
TopologyResidue* TopologyAssembly::GetResidueByIndex(int index)
{
    for(TopologyResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        TopologyResidue* residue = (*it);
        if(residue->GetIndex() == index)
            return residue;
    }
    return NULL;
}
int TopologyAssembly::GetAtomIndexByName(string atom_name)
{
    for(TopologyResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        TopologyResidue* residue = (*it);
        TopologyResidue::TopologyAtomVector atoms = residue->GetAtoms();
        for(TopologyResidue::TopologyAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            TopologyAtom* atom = (*it1);
            string name = atom->GetAtomName();
            if(name.compare(atom_name) == 0)
            {
                return atom->GetIndex();
            }
        }
    }
    return -1;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAssembly::SetAssemblyName(const string assembly_name)
{
    assembly_name_ = assembly_name;
}
void TopologyAssembly::SetResidues(TopologyResidueVector residues)
{
    residues_.clear();
    for(TopologyResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        TopologyResidue* residue = (*it);
        residues_.push_back(residue);
    }
}
void TopologyAssembly::AddResidue(TopologyResidue *residue)
{
    residues_.push_back(residue);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAssembly::Print(ostream &out)
{
    out << "Assembly Name: " << assembly_name_ << endl;
    for(TopologyResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        TopologyResidue* residue = (*it);
        out << "-------------------------- " << residue->GetResidueName() << " -------------------------" << endl;
        residue->Print(out);
    }
}

