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
#include "../../../../includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
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
PdbFile* Assembly::BuildPdbFileStructureFromAssembly(int link_card_direction, int connect_card_existance)
{
    cout << "Creating PDB file" << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating PDB file ...");
    PdbFile* pdb_file = new PdbFile();
    PdbTitleSection* title_card = new PdbTitleSection();
    title_card->SetTitle("Generated by GMML");
    // Set pdb_file title card
    pdb_file->SetTitle(title_card);

    PdbModelSection* model_card = new PdbModelSection();
    PdbModelSection::PdbModelCardMap models = PdbModelSection::PdbModelCardMap();
    PdbModelCard* model = new PdbModelCard();
    model->SetModelSerialNumber(1);
    PdbModelResidueSet* residue_set = new PdbModelResidueSet();
    int serial_number = 1;
    int sequence_number = 1;

    AssemblytoPdbSequenceNumberMap assembly_to_sequence_number_map = AssemblytoPdbSequenceNumberMap();
    AssemblytoPdbSerialNumberMap assembly_to_serial_number_map = AssemblytoPdbSerialNumberMap();
    ExtractPdbModelSectionFromAssembly(residue_set, serial_number, sequence_number, model_index_, assembly_to_sequence_number_map,
                                    assembly_to_serial_number_map);

    PdbLinkSection* link_card = new PdbLinkSection();
    ExtractPdbLinkSectionFromAssembly(link_card, model_index_, assembly_to_sequence_number_map, link_card_direction);
    link_card->SetRecordName("LINK");
    pdb_file->SetLinks(link_card);

    if(connect_card_existance == 1)
    {
        PdbConnectSection* connect_card = new PdbConnectSection();
        ExtractPdbConnectSectionFromAssembly(connect_card, assembly_to_serial_number_map);
        pdb_file->SetConnectivities(connect_card);
    }
    model->SetModelResidueSet(residue_set);
    models[1] = model;
    model_card->SetModels(models);
    pdb_file->SetModels(model_card);

    cout << "PDB file created" << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "PDB file created");
    return pdb_file;
}

void Assembly::ExtractPdbModelSectionFromAssembly(PdbModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number,
                                               AssemblytoPdbSequenceNumberMap& assembly_to_pdb_sequence_number_map,
                                               AssemblytoPdbSerialNumberMap& assembly_to_pdb_serial_number_map)
{
    for(AssemblyVector::iterator it = this->assemblies_.begin(); it != this->assemblies_.end(); it++)
    {
        Assembly* assembly = (*it);
        AssemblyVector assemblies = assembly->GetAssemblies();
        for(AssemblyVector::iterator it1 = assemblies.begin(); it1 != assemblies.end(); it1++)
        {
            ExtractPdbModelSectionFromAssembly(residue_set, serial_number, sequence_number, model_number, assembly_to_pdb_sequence_number_map,
                                            assembly_to_pdb_serial_number_map);
        }
        PdbAtomSection* atom_card = new PdbAtomSection();
        PdbHeterogenAtomSection* het_atom_card = new PdbHeterogenAtomSection();
        PdbAtomSection::PdbAtomCardMap atom_map = PdbAtomSection::PdbAtomCardMap();
        PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbAtomSection::PdbAtomCardOrderVector();
        PdbHeterogenAtomSection::PdbHeterogenAtomCardMap het_atom_map = PdbHeterogenAtomSection::PdbHeterogenAtomCardMap();
        PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector het_atom_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
        ResidueVector residues = assembly->GetResidues();
        for(ResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            Residue* residue = (*it1);
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                Atom* atom = (*it2);
                vector<string> dscr = Split(atom->GetDescription(), ";");
                //                stringstream ss;
                //                ss << setw(2) << fixed << setprecision(1) << atom->MolecularDynamicAtom::GetCharge();
                string residue_name = atom->GetResidue()->GetName();
                if(residue_name.compare("TIP3PBOX") == 0 || residue_name.compare("TIP5PBOX") == 0)
                    residue_name = "HOH";
                PdbAtomCard* pdb_atom = new PdbAtomCard(serial_number, atom->GetName(), ' ', residue_name, ' ', sequence_number, ' ',
                                                *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->GetElementSymbol(), "");//ss.str());
                vector<string> atom_id_tokens = Split(atom->GetId(), "_");
                pdb_atom->SetAtomChainId(ConvertString<char>(atom_id_tokens.at(3)));
                pdb_atom->SetAtomInsertionCode(ConvertString<char>(atom_id_tokens.at(5)));
                pdb_atom->SetAtomAlternateLocation(ConvertString<char>(atom_id_tokens.at(6)));
                assembly_to_pdb_sequence_number_map[gmml::ConvertString<int>(atom_id_tokens.at(4))] = sequence_number;
                assembly_to_pdb_serial_number_map[ConvertString<int>(atom_id_tokens.at(1))] = serial_number;

                if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
                {
                    atom_map[serial_number] = pdb_atom;
                    atom_vector.push_back(pdb_atom);
                    serial_number++;
                }
                else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
                {
                    het_atom_map[serial_number] = pdb_atom;
                    het_atom_vector.push_back(pdb_atom);
                    serial_number++;
                }
                else
                {
                    atom_map[serial_number] = pdb_atom;
                    atom_vector.push_back(pdb_atom);
                    serial_number++;
                }
            }
            sequence_number++;
        }
        atom_card->SetAtomCards(atom_map);
        atom_card->SetOrderedAtomCards(atom_vector);
        het_atom_card->SetHeterogenAtoms(het_atom_map);
        het_atom_card->SetOrderedHeterogenAtoms(het_atom_vector);
        residue_set->AddAtom(atom_card);
        residue_set->AddHeterogenAtom(het_atom_card);
    }
    PdbAtomSection* atom_card = new PdbAtomSection();
    PdbHeterogenAtomSection* het_atom_card = new PdbHeterogenAtomSection();
    PdbAtomSection::PdbAtomCardMap atom_map = PdbAtomSection::PdbAtomCardMap();
    PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbAtomSection::PdbAtomCardOrderVector();
    PdbHeterogenAtomSection::PdbHeterogenAtomCardMap het_atom_map = PdbHeterogenAtomSection::PdbHeterogenAtomCardMap();
    PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector het_atom_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
    for(ResidueVector::iterator it1 = residues_.begin(); it1 != residues_.end(); it1++)
    {
        Residue* residue = (*it1);
        AtomVector atoms = residue->GetAtoms();
        for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            Atom* atom = (*it2);
            vector<string> dscr = Split(atom->GetDescription(), ";");
            //            stringstream ss;
            //            ss << setw(2) << fixed << setprecision(1) << atom->MolecularDynamicAtom::GetCharge();
            string residue_name = atom->GetResidue()->GetName();
            if(residue_name.compare("TIP3PBOX") == 0 || residue_name.compare("TIP5PBOX") == 0)
                residue_name = "HOH";
            PdbAtomCard* pdb_atom = new PdbAtomCard(serial_number, atom->GetName(), ' ', residue_name, ' ', sequence_number, ' ',
                                            *((atom->GetCoordinates()).at(model_number)), dNotSet, dNotSet, atom->GetElementSymbol(), "");//ss.str());
            vector<string> atom_id_tokens = Split(atom->GetId(), "_");
            pdb_atom->SetAtomChainId(ConvertString<char>(atom_id_tokens.at(3)));
            pdb_atom->SetAtomInsertionCode(ConvertString<char>(atom_id_tokens.at(5)));
            pdb_atom->SetAtomAlternateLocation(ConvertString<char>(atom_id_tokens.at(6)));
            assembly_to_pdb_sequence_number_map[gmml::ConvertString<int>(atom_id_tokens.at(4))] = sequence_number;
            assembly_to_pdb_serial_number_map[ConvertString<int>(atom_id_tokens.at(1))] = serial_number;

            if(find(dscr.begin(), dscr.end(), "Atom") != dscr.end())
            {

                atom_map[serial_number] = pdb_atom;
                atom_vector.push_back(pdb_atom);
                serial_number++;
            }
            else if(find(dscr.begin(), dscr.end(), "Het") != dscr.end())
            {
                het_atom_map[serial_number] = pdb_atom;
                het_atom_vector.push_back(pdb_atom);
                serial_number++;
            }
            else
            {
                atom_map[serial_number] = pdb_atom;
                atom_vector.push_back(pdb_atom);
                serial_number++;
            }
        }
        sequence_number++;
    }
    atom_card->SetAtomCards(atom_map);
    atom_card->SetOrderedAtomCards(atom_vector);
    het_atom_card->SetHeterogenAtoms(het_atom_map);
    het_atom_card->SetOrderedHeterogenAtoms(het_atom_vector);
    residue_set->AddAtom(atom_card);
    residue_set->AddHeterogenAtom(het_atom_card);
}

void Assembly::ExtractPdbLinkSectionFromAssembly(PdbLinkSection* link_card, int model_index, AssemblytoPdbSequenceNumberMap assembly_to_pdb_sequence_number,
                                              int link_card_direction)
{
    PdbLinkSection::LinkCardVector link_vector = PdbLinkSection::LinkCardVector();
    vector<string> visited_links = vector<string>();
    AtomVector all_atoms = GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = (*it);
        Residue* residue = atom->GetResidue();
        AtomNode* node = atom->GetNode();
        if(node != NULL)
        {
            AtomVector neighbors = node->GetNodeNeighbors();

            if((link_card_direction == 1 && (atom->GetName().find("C") != string::npos) &&
                (find(visited_links.begin(), visited_links.end(), atom->GetId()) == visited_links.end())) ||
                    (link_card_direction == -1 && (atom->GetName().find("O") != string::npos) &&
                     (find(visited_links.begin(), visited_links.end(), atom->GetId()) == visited_links.end())))
            {
                for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
                {
                    Atom* neighbor = (*it1);
                    Residue* neighbor_residue = neighbor->GetResidue();
                    if(find(visited_links.begin(), visited_links.end(), neighbor->GetId()) == visited_links.end())
                    {
                        if(residue->GetId().compare(neighbor_residue->GetId()) != 0)
                        {
                            visited_links.push_back(atom->GetId());
                            visited_links.push_back(neighbor->GetId());
                            PdbLinkCardResidue* link_residue1 = new PdbLinkCardResidue();
                            vector<string> atom_id_tokens = Split(atom->GetId(), "_");
                            link_residue1->SetAtomName(atom->GetName());
                            link_residue1->SetResidueChainId(ConvertString<char>(atom_id_tokens.at(3)));
                            link_residue1->SetResidueSequenceNumber(assembly_to_pdb_sequence_number[ConvertString<int>(atom_id_tokens.at(4))]);
                            link_residue1->SetResidueInsertionCode(ConvertString<char>(atom_id_tokens.at(5)));
                            link_residue1->SetAlternateLocationIndicator(ConvertString<char>(atom_id_tokens.at(6)));
                            link_residue1->SetResidueName(residue->GetName());
                            PdbLinkCardResidue* link_residue2 = new PdbLinkCardResidue();
                            vector<string> neighbor_id_tokens = Split(neighbor->GetId(), "_");
                            link_residue2->SetAtomName(neighbor->GetName());
                            link_residue2->SetResidueChainId(ConvertString<char>(neighbor_id_tokens.at(3)));
                            link_residue2->SetResidueSequenceNumber(assembly_to_pdb_sequence_number[ConvertString<int>(neighbor_id_tokens.at(4))]);
                            link_residue2->SetResidueInsertionCode(ConvertString<char>(neighbor_id_tokens.at(5)));
                            link_residue2->SetAlternateLocationIndicator(ConvertString<char>(neighbor_id_tokens.at(6)));
                            link_residue2->SetResidueName(neighbor_residue->GetName());

                            PdbLinkCard* pdb_link = new PdbLinkCard();
                            pdb_link->AddResidue(link_residue1);
                            pdb_link->AddResidue(link_residue2);
                            double distance = atom->GetCoordinates().at(model_index)->Distance(*(neighbor->GetCoordinates().at(model_index)));
                            pdb_link->SetLinkLength(distance);

                            link_vector.push_back(pdb_link);
                        }
                    }
                }
            }
        }
    }
    link_card->SetResidueLinkCards(link_vector);
}

void Assembly::ExtractPdbConnectSectionFromAssembly(PdbConnectSection *connect_card, AssemblytoPdbSerialNumberMap assembly_to_pdb_serial_number)
{
    PdbConnectSection::BondedAtomsSerialNumbersMap bonded_atoms_serial_number_map = PdbConnectSection::BondedAtomsSerialNumbersMap();
    AtomVector all_atoms = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = *it;
        AtomNode* node = atom->GetNode();
        if(node != NULL)
        {
            AtomVector neighbors = node->GetNodeNeighbors();
            vector<string> atom_id_tokens = Split(atom->GetId(), "_");
            int atom_serial_number = assembly_to_pdb_serial_number[ConvertString<int>(atom_id_tokens.at(1))];
            bonded_atoms_serial_number_map[atom_serial_number] = vector<int>();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = *it1;
                vector<string> neighbor_id_tokens = Split(neighbor->GetId(), "_");
                bonded_atoms_serial_number_map[atom_serial_number].push_back(assembly_to_pdb_serial_number[ConvertString<int>(neighbor_id_tokens.at(1))]);
            }
        }
    }
    connect_card->SetBondedAtomsSerialNumbers(bonded_atoms_serial_number_map);
}
