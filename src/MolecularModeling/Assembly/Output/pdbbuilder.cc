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

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
PdbFileSpace::PdbFile* Assembly::BuildPdbFileStructureFromAssembly(int link_card_direction, int connect_card_existance, int model_index)
{
    if (model_index == -1) // -1 is the default set in header file, if that's not changed, then use model_index_. Can't pass in model_index_ as default for reasons.
        model_index = model_index_;

//    std::cout << "Creating PDB file" << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating PDB file ...");
    PdbFileSpace::PdbFile* pdb_file = new PdbFileSpace::PdbFile();
    PdbFileSpace::PdbTitleSection* title_card = new PdbFileSpace::PdbTitleSection();
    title_card->SetTitle("Generated by GMML");
    // Set pdb_file title card
    pdb_file->SetTitle(title_card);

    PdbFileSpace::PdbModelSection* model_card = new PdbFileSpace::PdbModelSection();
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = PdbFileSpace::PdbModelSection::PdbModelCardMap();
    PdbFileSpace::PdbModelCard* model = new PdbFileSpace::PdbModelCard();
    model->SetModelSerialNumber(1);
    PdbFileSpace::PdbModelResidueSet* residue_set = new PdbFileSpace::PdbModelResidueSet();
    int serial_number = 1;
    int sequence_number = 1;

    AssemblytoPdbSequenceNumberMap assembly_to_sequence_number_map = AssemblytoPdbSequenceNumberMap();
    AssemblytoPdbSerialNumberMap assembly_to_serial_number_map = AssemblytoPdbSerialNumberMap();
    ExtractPdbModelSectionFromAssembly(residue_set, serial_number, sequence_number, model_index, assembly_to_sequence_number_map,
                                    assembly_to_serial_number_map);

    PdbFileSpace::PdbLinkSection* link_card = new PdbFileSpace::PdbLinkSection();
    //The follwing line might be commented out for my testing purpose, definitely shouldn'be committed/pushed. If you see it commented out, please uncomment it.
    ExtractPdbLinkSectionFromAssembly(link_card, model_index, assembly_to_sequence_number_map, link_card_direction);
    link_card->SetRecordName("LINK");
    pdb_file->SetLinks(link_card);

    if(connect_card_existance == 1)  //The CONECT card generated by this code cannot be read correctly by Chimera. Cannot find out the bug. --Yao
    {
        PdbFileSpace::PdbConnectSection* connect_card = new PdbFileSpace::PdbConnectSection();
        ExtractPdbConnectSectionFromAssembly(connect_card, assembly_to_serial_number_map);
        pdb_file->SetConnectivities(connect_card);
    }
    model->SetModelResidueSet(residue_set);
    models[1] = model;
    model_card->SetModels(models);
    pdb_file->SetModels(model_card);

    gmml::log(__LINE__, __FILE__, gmml::INF, "PDB file created");
    return pdb_file;
}

void Assembly::ExtractPdbModelSectionFromAssembly(PdbFileSpace::PdbModelResidueSet* residue_set, int &serial_number, int &sequence_number, int model_number,
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
        PdbFileSpace::PdbAtomSection* atom_card = new PdbFileSpace::PdbAtomSection();
        PdbFileSpace::PdbHeterogenAtomSection* het_atom_card = new PdbFileSpace::PdbHeterogenAtomSection();
        PdbFileSpace::PdbAtomSection::PdbAtomMap atom_map = PdbFileSpace::PdbAtomSection::PdbAtomMap();
        PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector();
        PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomCardMap het_atom_map = PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomCardMap();
        PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector het_atom_vector = PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
        //ResidueVector residues = assembly->GetResidues();
        ResidueVector residues = this->GetResidues();
        for(ResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            Residue* residue = (*it1);
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                Atom* atom = (*it2);
                std::vector<std::string> dscr = gmml::Split(atom->GetDescription(), ";");
                //                std::stringstream ss;
                //                ss << setw(2) << fixed << setprecision(1) << atom->MolecularDynamicAtom::GetCharge();
                std::string residue_name = atom->GetResidue()->GetName();
                if(residue_name.compare("TIP3PBOX") == 0 || residue_name.compare("TIP5PBOX") == 0)
                    residue_name = "HOH";
                PdbFileSpace::PdbAtomCard* pdb_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetName(), ' ', residue_name, ' ', sequence_number, ' ',
                                                *((atom->GetCoordinates()).at(model_number)), gmml::dNotSet, gmml::dNotSet, atom->GetElementSymbol(), "");//ss.str());
                std::vector<std::string> atom_id_tokens = gmml::Split(atom->GetId(), "_");
                pdb_atom->SetAtomChainId(gmml::ConvertString<char>(atom_id_tokens.at(3)));
                pdb_atom->SetAtomInsertionCode(gmml::ConvertString<char>(atom_id_tokens.at(5)));
                pdb_atom->SetAtomAlternateLocation(gmml::ConvertString<char>(atom_id_tokens.at(6)));
                assembly_to_pdb_sequence_number_map[gmml::ConvertString<int>(atom_id_tokens.at(4))] = sequence_number;
                assembly_to_pdb_serial_number_map[gmml::ConvertString<int>(atom_id_tokens.at(1))] = serial_number;

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
    PdbFileSpace::PdbAtomSection* atom_card = new PdbFileSpace::PdbAtomSection();
    PdbFileSpace::PdbHeterogenAtomSection* het_atom_card = new PdbFileSpace::PdbHeterogenAtomSection();
    PdbFileSpace::PdbAtomSection::PdbAtomMap atom_map = PdbFileSpace::PdbAtomSection::PdbAtomMap();
    PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector();
    PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomCardMap het_atom_map = PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomCardMap();
    PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector het_atom_vector = PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
    for(ResidueVector::iterator it1 = residues_.begin(); it1 != residues_.end(); it1++)
    {
        Residue* residue = (*it1);
        AtomVector atoms = residue->GetAtoms();
        for(AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            Atom* atom = (*it2);
            std::vector<std::string> dscr = gmml::Split(atom->GetDescription(), ";");
            //            std::stringstream ss;
            //            ss << setw(2) << fixed << setprecision(1) << atom->MolecularDynamicAtom::GetCharge();
            std::string residue_name = atom->GetResidue()->GetName();
            if(residue_name.compare("TIP3PBOX") == 0 || residue_name.compare("TIP5PBOX") == 0)
                residue_name = "HOH";
            PdbFileSpace::PdbAtomCard* pdb_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetName(), ' ', residue_name, ' ', sequence_number, ' ',
                                            *((atom->GetCoordinates()).at(model_number)), gmml::dNotSet, gmml::dNotSet, atom->GetElementSymbol(), "");//ss.str());
            std::vector<std::string> atom_id_tokens = gmml::Split(atom->GetId(), "_");
            pdb_atom->SetAtomChainId(gmml::ConvertString<char>(atom_id_tokens.at(3)));
            pdb_atom->SetAtomInsertionCode(gmml::ConvertString<char>(atom_id_tokens.at(5)));
            pdb_atom->SetAtomAlternateLocation(gmml::ConvertString<char>(atom_id_tokens.at(6)));
            assembly_to_pdb_sequence_number_map[gmml::ConvertString<int>(atom_id_tokens.at(4))] = sequence_number;
            assembly_to_pdb_serial_number_map[gmml::ConvertString<int>(atom_id_tokens.at(1))] = serial_number;

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

void Assembly::ExtractPdbLinkSectionFromAssembly(PdbFileSpace::PdbLinkSection* link_card, int model_index, AssemblytoPdbSequenceNumberMap assembly_to_pdb_sequence_number,
                                              int link_card_direction)
{
    PdbFileSpace::PdbLinkSection::LinkCardVector link_vector = PdbFileSpace::PdbLinkSection::LinkCardVector();
    std::vector<std::string> visited_links = std::vector<std::string>();
    AtomVector all_atoms = GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = (*it);
        Residue* residue = atom->GetResidue();
        AtomNode* node = atom->GetNode();
        if(node != NULL)
        {
            AtomVector neighbors = node->GetNodeNeighbors();

            if((link_card_direction == 1 && (atom->GetName().find("C") != std::string::npos) &&
                (find(visited_links.begin(), visited_links.end(), atom->GetId()) == visited_links.end())) ||
                    (link_card_direction == -1 && (atom->GetName().find("O") != std::string::npos) &&
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
			    //my code
				if (residue->GetId().empty()) {
				}
			    std::vector <std::string> ResidueIdSplit = gmml::Split(residue->GetId(),"_");
			    std::vector <std::string> NeighborResidueIdSplit = gmml::Split(neighbor_residue->GetId(),"_");
			    std::string ResidueNumber = ResidueIdSplit.at(2);
			    std::string NeighborResidueNumber = NeighborResidueIdSplit.at(2);
			
			    if (ResidueNumber.compare(NeighborResidueNumber) != 0) // if this pair is not the two partial occupancies of the same residue
			    {//my code
                            visited_links.push_back(atom->GetId());
                            visited_links.push_back(neighbor->GetId());
                            PdbFileSpace::PdbLinkCardResidue* link_residue1 = new PdbFileSpace::PdbLinkCardResidue();
                            std::vector<std::string> atom_id_tokens = gmml::Split(atom->GetId(), "_");
                            link_residue1->SetAtomName(atom->GetName());
                            link_residue1->SetResidueChainId(gmml::ConvertString<char>(atom_id_tokens.at(3)));
                            link_residue1->SetResidueSequenceNumber(assembly_to_pdb_sequence_number[gmml::ConvertString<int>(atom_id_tokens.at(4))]);
                            link_residue1->SetResidueInsertionCode(gmml::ConvertString<char>(atom_id_tokens.at(5)));
                            link_residue1->SetAlternateLocationIndicator(gmml::ConvertString<char>(atom_id_tokens.at(6)));
                            link_residue1->SetResidueName(residue->GetName());
                            PdbFileSpace::PdbLinkCardResidue* link_residue2 = new PdbFileSpace::PdbLinkCardResidue();
                            std::vector<std::string> neighbor_id_tokens = gmml::Split(neighbor->GetId(), "_");
                            link_residue2->SetAtomName(neighbor->GetName());
                            link_residue2->SetResidueChainId(gmml::ConvertString<char>(neighbor_id_tokens.at(3)));
                            link_residue2->SetResidueSequenceNumber(assembly_to_pdb_sequence_number[gmml::ConvertString<int>(neighbor_id_tokens.at(4))]);
                            link_residue2->SetResidueInsertionCode(gmml::ConvertString<char>(neighbor_id_tokens.at(5)));
                            link_residue2->SetAlternateLocationIndicator(gmml::ConvertString<char>(neighbor_id_tokens.at(6)));
                            link_residue2->SetResidueName(neighbor_residue->GetName());

                            PdbFileSpace::PdbLinkCard* pdb_link = new PdbFileSpace::PdbLinkCard();
                            pdb_link->AddResidue(link_residue1);
                            pdb_link->AddResidue(link_residue2);
                            double distance = atom->GetCoordinates().at(model_index)->Distance(*(neighbor->GetCoordinates().at(model_index)));
                            pdb_link->SetLinkLength(distance);

                            link_vector.push_back(pdb_link);
			    }//my code
                        }
                    }
                }
            }
        }
    }
    link_card->SetResidueLinkCards(link_vector);
}

void Assembly::ExtractPdbConnectSectionFromAssembly(PdbFileSpace::PdbConnectSection *connect_card, AssemblytoPdbSerialNumberMap assembly_to_pdb_serial_number)
{
    PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap bonded_atoms_serial_number_map = PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap();
    AtomVector all_atoms = this->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = all_atoms.begin(); it != all_atoms.end(); it++)
    {
        Atom* atom = *it;
        AtomNode* node = atom->GetNode();
        if(node != NULL)
        {
            AtomVector neighbors = node->GetNodeNeighbors();
            std::vector<std::string> atom_id_tokens = gmml::Split(atom->GetId(), "_");
            int atom_serial_number = assembly_to_pdb_serial_number[gmml::ConvertString<int>(atom_id_tokens.at(1))];
            bonded_atoms_serial_number_map[atom_serial_number] = std::vector<int>();
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = *it1;
                std::vector<std::string> neighbor_id_tokens = gmml::Split(neighbor->GetId(), "_");
                bonded_atoms_serial_number_map[atom_serial_number].push_back(assembly_to_pdb_serial_number[gmml::ConvertString<int>(neighbor_id_tokens.at(1))]);
            }
        }
    }
    connect_card->SetBondedAtomsSerialNumbers(bonded_atoms_serial_number_map);
}
