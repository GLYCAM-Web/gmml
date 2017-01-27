#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblink.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/grid.hpp"
#include "../../../includes/GeometryTopology/cell.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>

using namespace std;
using namespace MolecularModeling;
using namespace TopologyFileSpace;
using namespace CoordinateFileSpace;
using namespace PrepFileSpace;
using namespace PdbFileSpace;
using namespace PdbqtFileSpace;
using namespace ParameterFileSpace;
using namespace GeometryTopology;
using namespace LibraryFileSpace;
using namespace gmml;
using namespace Glycan;
using namespace CondensedSequenceSpace;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Assembly::AttachResidues(Residue *residue, Residue *parent_residue, int branch_index, string parameter_file)
{
    ///Translate all atoms of the attached residue to place them in proper position with respect to the tail atom of the parent residue/assembly
    this->SetAttachedResidueBond(residue, parent_residue, branch_index, parameter_file);

    ///Rotate all atoms of the attached residue to set the proper bond angle between the attached residue and the parent residue
    this->SetAttachedResidueAngle(residue, parent_residue, branch_index, parameter_file);

    ///Rotate all atoms of the attached residue to set the proper Phi, Psi and Omega torsion angles
    this->SetAttachedResidueTorsion(residue, parent_residue, branch_index);
}

void Assembly::RemoveHydrogenAtAttachedPosition(Residue *residue, int branch_index)
{
    Atom* oxygen = residue->GetTailAtoms().at(branch_index);
    if(oxygen != NULL)
    {
        int oxygen_index = 1;
        if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
            oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));

        Atom* hydrogen = NULL;

        AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
        {
            Atom* neighbor = *it;
            if(neighbor->GetName().at(0) == 'H' &&
                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
            {
                hydrogen = neighbor;
                break;
            }

        }
        residue->RemoveAtom(hydrogen);
    }
}

void Assembly::AdjustCharge(Residue *residue, Residue *parent_residue, int branch_index)
{
    if(residue->GetName().compare("SO3") == 0)
    {
      parent_residue->GetTailAtoms().at(branch_index)->MolecularDynamicAtom::SetCharge(
                  parent_residue->GetTailAtoms().at(branch_index)->MolecularDynamicAtom::GetCharge() + 0.031);
    }
    else if(residue->GetName().compare("MEX") == 0 || residue->GetName().compare("ACX") == 0)
    {
        Atom* oxygen = parent_residue->GetTailAtoms().at(branch_index);
        Atom* carbon = NULL;
        int oxygen_index = 1;
        if(oxygen->GetName().size() > 1 && isdigit(oxygen->GetName().at(1)))
            oxygen_index = ConvertString<int>(ConvertT<char>(oxygen->GetName().at(1)));

        AtomVector oxygen_neighbors = oxygen->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = oxygen_neighbors.begin(); it != oxygen_neighbors.end(); it++)
        {
            Atom* neighbor = *it;
            if(neighbor->GetName().at(0) == 'C' &&
                    (neighbor->GetName().size() > 1 && isdigit(neighbor->GetName().at(1)) &&
                     ConvertString<int>(ConvertT<char>(neighbor->GetName().at(1))) == oxygen_index))
            {
                carbon = neighbor;
                break;
            }
        }
        if(carbon != NULL)
        {
            if(residue->GetName().compare("MEX") == 0)
                carbon->MolecularDynamicAtom::SetCharge(carbon->MolecularDynamicAtom::GetCharge() - 0.039);
            if(residue->GetName().compare("ACX") == 0)
                carbon->MolecularDynamicAtom::SetCharge(carbon->MolecularDynamicAtom::GetCharge() + 0.008);
        }
    }
}
