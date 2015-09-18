#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModelResidueSet::PdbqtModelResidueSet()
{
    roots_ = NULL;
    atoms_ = NULL;
}

PdbqtModelResidueSet::PdbqtModelResidueSet(stringstream &residue_set_block)
{
    roots_ = NULL;
    atoms_ = NULL;
    branches_ = BranchCardVector();
    string line;
    getline(residue_set_block, line);
    string temp = line;
    stringstream root_block;
    stringstream all_atoms_block;
    if(line.find("ROOT") != string::npos)
    {
        root_block << line << endl;
        getline(residue_set_block, line);
        while(line.find("ENDROOT") == string::npos)
        {
            if(line.find("ATOM") != string::npos || line.find("HETATM") != string::npos)
            {
                root_block << line << endl;
                all_atoms_block << line << endl;
            }
            getline(residue_set_block, line);
        }
        root_block << line << endl;
    }
    roots_ = new PdbqtRootCard(root_block);
    getline(residue_set_block, line);
    while(!Trim(temp).empty())
    {
        stringstream branch_block;
        if(line.find("BRANCH") != string::npos)
        {
            branch_block << line << endl;
            vector<string> tokens = Split(line, " ");
            int solid_atom_serial_number = ConvertString<int>(tokens.at(1));
            int rotatable_atom_serial_number = ConvertString<int>(tokens.at(2));
            stringstream end_branch_line;
            end_branch_line << "ENDBRANCH " << setw(3) << solid_atom_serial_number << " " << setw(3) << rotatable_atom_serial_number;
            getline(residue_set_block, line);
            while(line.find(end_branch_line.str()) == string::npos)
            {
                if(line.find("ATOM") != string::npos || line.find("HETATM") != string::npos)
                {
                    all_atoms_block << line << endl;
                }
                branch_block << line << endl;
                getline(residue_set_block, line);
            }
            branch_block << line << endl;
        }
        branches_.push_back(new PdbqtBranchCard(branch_block));
        getline(residue_set_block, line);
        if(line.empty())
            break;
    }
    atoms_ = new PdbqtAtomCard(all_atoms_block);
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
PdbqtRootCard* PdbqtModelResidueSet::GetRoots()
{
    return roots_;
}
PdbqtModelResidueSet::BranchCardVector PdbqtModelResidueSet::GetBranches()
{
    return branches_;
}
PdbqtAtomCard* PdbqtModelResidueSet::GetAtoms()
{
    return atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtModelResidueSet::SetRoots(PdbqtRootCard* roots)
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
void PdbqtModelResidueSet::SetAtoms(PdbqtAtomCard *atoms)
{
    atoms_ = atoms;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtModelResidueSet::Print(ostream &out)
{
    out << "====================== Roots =====================" << endl;
    if(roots_ != NULL)
        roots_->Print(out);
    for(BranchCardVector::iterator it = branches_.begin(); it != branches_.end(); it++)
    {
        out << "====================== Main Branch =====================" << endl;
        (*it)->Print(out);
    }
}


