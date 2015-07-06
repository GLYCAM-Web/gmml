#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbqtFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModelResidueSet::PdbqtModelResidueSet(){}

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

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtModelResidueSet::Print(ostream &out)
{
}


