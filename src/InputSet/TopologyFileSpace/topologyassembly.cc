
#include "../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"

using TopologyFileSpace::TopologyAssembly;

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
std::string TopologyAssembly::GetAssemblyName()
{
    return assembly_name_;
}
TopologyAssembly::TopologyResidueVector TopologyAssembly::GetResidues()
{
    return residues_;
}
TopologyFileSpace::TopologyResidue* TopologyAssembly::GetResidueByIndex(int index)
{
    for(TopologyResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        TopologyFileSpace::TopologyResidue* residue = (*it);
        if(residue->GetIndex() == index)
            return residue;
    }
    return NULL;
}
int TopologyAssembly::GetAtomIndexByName(std::string atom_name, int residue_index)
{
    for(TopologyResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        TopologyFileSpace::TopologyResidue* residue = (*it);
        if(residue->GetIndex() == residue_index)
        {
            TopologyResidue::TopologyAtomVector atoms = residue->GetAtoms();
            for(TopologyFileSpace::TopologyResidue::TopologyAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
            {
                TopologyAtom* atom = (*it1);
                std::string name = atom->GetAtomName();
                if(name.compare(atom_name) == 0)
                {
                    return atom->GetIndex();
                }
            }
        }
    }
    return -1;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void TopologyAssembly::SetAssemblyName(const std::string assembly_name)
{
    assembly_name_ = assembly_name;
}
void TopologyAssembly::SetResidues(TopologyResidueVector residues)
{
    residues_.clear();
    for(TopologyResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        TopologyFileSpace::TopologyResidue* residue = (*it);
        residues_.push_back(residue);
    }
}
void TopologyAssembly::AddResidue(TopologyFileSpace::TopologyResidue *residue)
{
    residues_.push_back(residue);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void TopologyAssembly::Print(std::ostream &out)
{
    out << "Assembly Name: " << assembly_name_ << std::endl;
    for(TopologyResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        TopologyFileSpace::TopologyResidue* residue = (*it);
        out << "-------------------------- " << residue->GetResidueName() << " -------------------------" << std::endl;
        residue->Print(out);
    }
}
