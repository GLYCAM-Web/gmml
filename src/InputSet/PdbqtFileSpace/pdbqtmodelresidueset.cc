#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtModelResidueSet;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModelResidueSet::PdbqtModelResidueSet()
{
    roots_ = NULL;
    atoms_ = NULL;
}

PdbqtModelResidueSet::PdbqtModelResidueSet(std::stringstream &residue_set_block)
{
    roots_ = NULL;
    atoms_ = NULL;
    branches_ = BranchCardVector();
    std::string line;
    getline(residue_set_block, line);
    std::string temp = line;
    std::stringstream root_block;
    std::stringstream all_atoms_block;
    if(line.find("ROOT") != std::string::npos)
    {
        root_block << line << std::endl;
        getline(residue_set_block, line);
        while(line.find("ENDROOT") == std::string::npos)
        {
            if(line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos)
            {
                root_block << line << std::endl;
                all_atoms_block << line << std::endl;
            }
            getline(residue_set_block, line);
        }
        root_block << line << std::endl;
    }
    roots_ = new PdbqtFileSpace::PdbqtRootCard(root_block);
    getline(residue_set_block, line);
    while(!gmml::Trim(temp).empty())
    {
        std::stringstream branch_block;
        if(line.find("BRANCH") != std::string::npos)
        {
            branch_block << line << std::endl;
            std::vector<std::string> tokens = gmml::Split(line, " ");
            int solid_atom_serial_number = gmml::ConvertString<int>(tokens.at(1));
            int rotatable_atom_serial_number = gmml::ConvertString<int>(tokens.at(2));
            std::stringstream end_branch_line;
            end_branch_line << "ENDBRANCH " << std::setw(3) << solid_atom_serial_number << " " << std::setw(3) << rotatable_atom_serial_number;
            getline(residue_set_block, line);
            while(line.find(end_branch_line.str()) == std::string::npos)
            {
                if(line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos)
                {
                    all_atoms_block << line << std::endl;
                }
                branch_block << line << std::endl;
                getline(residue_set_block, line);
            }
            branch_block << line << std::endl;
        }
        branches_.push_back(new PdbqtBranchCard(branch_block));
        getline(residue_set_block, line);
        if(line.empty())
            break;
    }
    atoms_ = new PdbqtFileSpace::PdbqtAtomCard(all_atoms_block);
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
PdbqtFileSpace::PdbqtRootCard* PdbqtModelResidueSet::GetRoots()
{
    return roots_;
}
PdbqtModelResidueSet::BranchCardVector PdbqtModelResidueSet::GetBranches()
{
    return branches_;
}
PdbqtFileSpace::PdbqtAtomCard* PdbqtModelResidueSet::GetAtoms()
{
    return atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtModelResidueSet::SetRoots(PdbqtFileSpace::PdbqtRootCard* roots)
{
    roots_ = roots;
}
void PdbqtModelResidueSet::SetBranches(BranchCardVector branches)
{
    branches_.clear();
    for(BranchCardVector::iterator it = branches.begin(); it != branches.end(); it++)
    {
        branches_.push_back(*it);
    }
}
void PdbqtModelResidueSet::AddBranch(PdbqtBranchCard* branch)
{
    branches_.push_back(branch);
}
void PdbqtModelResidueSet::SetAtoms(PdbqtFileSpace::PdbqtAtomCard *atoms)
{
    atoms_ = atoms;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtModelResidueSet::Print(std::ostream &out)
{
    out << "====================== Roots =====================" << std::endl;
    if(roots_ != NULL)
        roots_->Print(out);
    for(BranchCardVector::iterator it = branches_.begin(); it != branches_.end(); it++)
    {
        out << "====================== Main Branch =====================" << std::endl;
        (*it)->Print(out);
    }
}
