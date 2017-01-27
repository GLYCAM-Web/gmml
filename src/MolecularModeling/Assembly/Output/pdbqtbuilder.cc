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
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
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
#include "../../../../includes/InputSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblink.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
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
PdbqtFile* Assembly::BuildPdbqtFileStructureFromAssembly()
{
    cout << "Creating PDBQT file" << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating PDBQT file ...");
    PdbqtFile* pdbqt_file = new PdbqtFile();

    PdbqtModelCard* model_card = new PdbqtModelCard();
    PdbqtModelCard::PdbqtModelMap models = PdbqtModelCard::PdbqtModelMap();
    PdbqtModel* model = new PdbqtModel();
    model->SetModelSerialNumber(1);
    PdbqtModelResidueSet* residue_set = new PdbqtModelResidueSet();
    int serial_number = 1;
    int sequence_number = 1;

    ExtractPdbqtModelCardFromAssembly(residue_set, serial_number, sequence_number, model_index_);

    model->SetModelResidueSet(residue_set);
    models[1] = model;
    model_card->SetModels(models);
    pdbqt_file->SetModels(model_card);
    cout << "PDBQT file created" << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "PDBQT file created");

    return pdbqt_file;
}

void Assembly::ExtractPdbqtModelCardFromAssembly(PdbqtModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number)
{
    PdbqtAtomCard* atom_card = new PdbqtAtomCard();
    PdbqtAtomCard::PdbqtAtomMap atom_map = PdbqtAtomCard::PdbqtAtomMap();
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        AssemblyVector assemblies = assembly->GetAssemblies();
        for(AssemblyVector::iterator it1 = assemblies.begin(); it1 != assemblies.end(); it1++)
        {
            ExtractPdbqtModelCardFromAssembly(residue_set, serial_number, sequence_number, model_number);
        }
        ResidueVector residues = assembly->GetResidues();
        for(ResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            Residue* residue = (*it1);
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                Atom* atom = (*it2);
                vector<string> dscr = Split(atom->GetDescription(), ";");
                if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
                {
                    PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                        *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                        atom->GetAtomType(), "ATOM");

                    atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
                else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
                {
                    PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                        *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                        atom->GetAtomType(), "HETATOM");
                    atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
                else
                {
                    PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                        *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                        atom->GetAtomType(), "ATOM");

                    atom_map[serial_number] = pdb_atom;
                    serial_number++;
                }
            }
            sequence_number++;
        }
    }
    for(ResidueVector::iterator it1 = residues_.begin(); it1 != residues_.end(); it1++)
    {
        Residue* residue = (*it1);
        AtomVector atoms = residue->GetAtoms();
        for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            Atom* atom = (*it2);
            vector<string> dscr = Split(atom->GetDescription(), ";");

            if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
            {
                PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                    *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                    atom->GetAtomType(), "ATOM");
                atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
            else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
            {
                PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                    *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                    atom->GetAtomType(), "HETATOM");
                atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
            else
            {
                PdbqtAtom* pdb_atom = new PdbqtAtom(serial_number, atom->GetName(), ' ', atom->GetResidue()->GetName(), ' ', sequence_number, ' ',
                                                    *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->MolecularDynamicAtom::GetCharge(),
                                                    atom->GetAtomType(), "ATOM");
                atom_map[serial_number] = pdb_atom;
                serial_number++;
            }
        }
        sequence_number++;
    }
    atom_card->SetAtoms(atom_map);
    residue_set->SetAtoms(atom_card);
}

