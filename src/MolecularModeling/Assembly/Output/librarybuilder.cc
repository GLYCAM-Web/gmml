#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../../../includes/MolecularModeling/assembly.hpp"
#include "../../../../includes/MolecularModeling/residue.hpp"
#include "../../../../includes/MolecularModeling/atom.hpp"
#include "../../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/utils.hpp"
#include "../../../../includes/common.hpp"
#include "../../../../includes/GeometryTopology/grid.hpp"
#include "../../../../includes/GeometryTopology/cell.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
LibraryFileSpace::LibraryFile* Assembly::BuildLibraryFileStructureFromAssembly()
{
//    std::cout << "Creating library file ..." << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating library file ...");
    LibraryFileSpace::LibraryFile* library_file = new LibraryFileSpace::LibraryFile();
    LibraryFileSpace::LibraryFile::ResidueMap residue_map = LibraryFileSpace::LibraryFile::ResidueMap();
    ResidueVector residues_of_assembly = this->GetAllResiduesOfAssembly();
    for(ResidueVector::iterator it = residues_of_assembly.begin(); it != residues_of_assembly.end(); it++)
    {
        Residue* assembly_residue = *it;
        //        std::cout << assembly_residue->GetId() << std::endl;
        int residue_index = distance(residues_of_assembly.begin(), it) + 1;
        AtomVector assembly_residue_atoms = assembly_residue->GetAtoms();
        LibraryFileSpace::LibraryFileResidue* library_residue = new LibraryFileSpace::LibraryFileResidue();
        library_residue->SetName(assembly_residue->GetName());
        AtomVector head_atoms = assembly_residue->GetHeadAtoms();
        AtomVector tail_atoms = assembly_residue->GetTailAtoms();
        if(!head_atoms.empty())
        {
            int head_atom_index = distance(assembly_residue_atoms.begin(), find(assembly_residue_atoms.begin(), assembly_residue_atoms.end(), head_atoms.at(0))) + 1;
            library_residue->SetHeadAtomIndex(head_atom_index);
        }
        if(!tail_atoms.empty())
        {
            int tail_atom_index = distance(assembly_residue_atoms.begin(), find(assembly_residue_atoms.begin(), assembly_residue_atoms.end(), tail_atoms.at(0))) + 1;
            library_residue->SetTailAtomIndex(tail_atom_index);
        }
        int order = 1;
        for(AtomVector::iterator it1 = assembly_residue_atoms.begin(); it1 != assembly_residue_atoms.end(); it1++)
        {
            Atom* residue_atom = (*it1);
            int atom_index = distance(assembly_residue_atoms.begin(), it1) + 1;
            std::vector<int> bonded_atom_indices = std::vector<int>();
            AtomNode* atom_node = residue_atom->GetNode();
            if(atom_node != NULL)
            {
                AtomVector atom_neighbours = atom_node->GetNodeNeighbors();
                for(AtomVector::iterator it2 = atom_neighbours.begin(); it2 != atom_neighbours.end(); it2++)
                {
                    Atom* atom = *it2;
                    int bonded_atom_index = distance(assembly_residue_atoms.begin(), find(assembly_residue_atoms.begin(), assembly_residue_atoms.end(),
                                                                                          atom));
                    bonded_atom_indices.push_back(bonded_atom_index);
                }
            }
            LibraryFileSpace::LibraryFileAtom* atom = new LibraryFileSpace::LibraryFileAtom(residue_atom->GetAtomType(), residue_atom->GetName(), residue_index, atom_index,
                                                        gmml::iNotSet, residue_atom->MolecularDynamicAtom::GetCharge(),
                                                        *(residue_atom->GetCoordinates()[assembly_residue->GetAssembly()->GetModelIndex()]), bonded_atom_indices,
                    order);
            order++;
            library_residue->AddAtom(atom);
        }
        residue_map[library_residue->GetName()] = library_residue;
    }
    library_file->SetResidues(residue_map);
    return library_file;
}
