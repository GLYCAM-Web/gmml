#include <string>

#include "includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "includes/InputSet/PdbFileSpace/pdbresidue.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "includes/common.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp"

using PdbFileSpace::PdbFile;

void PdbFile::DeleteResidue(PdbResidue *residue)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms = PdbAtomSection::PdbAtomMap();
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back((atom));
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }

            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
    this->UpdateConnectCard();
}

void PdbFile::DeleteResidues(PdbResidueVector target_residues)
{
    std::map<std::string, PdbResidue* > pdb_residue_map = std::map<std::string, PdbResidue* >();
    for(PdbResidueVector::iterator it = target_residues.begin(); it != target_residues.end(); it++)
    {
        PdbResidue* residue = (*it);
        std::string target_residue_name = residue->GetResidueName();
        char target_residue_chain_id = residue->GetResidueChainId();
        int target_residue_sequence_number = residue->GetResidueSequenceNumber();
        char target_residue_insertion_code = residue->GetResidueInsertionCode();
        char target_residue_alternate_location = residue->GetResidueAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        std::string target_key = ss.str();
        pdb_residue_map[target_key] = residue;
    }

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms = PdbAtomSection::PdbAtomMap();
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_residue_map.find(key) == pdb_residue_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_residue_map.find(key) == pdb_residue_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
    this->UpdateConnectCard();
}

void PdbFile::DeleteResidueWithTheGivenModelNumber(PdbResidue *residue, int model_number)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms = PdbAtomSection::PdbAtomMap();
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;

        models_->SetModels(models);
        this->UpdateConnectCard();
    }
}

void PdbFile::DeleteResiduesWithTheGivenModelNumber(PdbResidueVector target_residues, int model_number)
{
    std::map<std::string, PdbResidue* > pdb_residue_map = std::map<std::string, PdbResidue* >();
    for(PdbResidueVector::iterator it = target_residues.begin(); it != target_residues.end(); it++)
    {
        PdbResidue* residue = (*it);
        std::string target_residue_name = residue->GetResidueName();
        char target_residue_chain_id = residue->GetResidueChainId();
        int target_residue_sequence_number = residue->GetResidueSequenceNumber();
        char target_residue_insertion_code = residue->GetResidueInsertionCode();
        char target_residue_alternate_location = residue->GetResidueAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        std::string target_key = ss.str();
        pdb_residue_map[target_key] = residue;
    }

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms = PdbAtomSection::PdbAtomMap();
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_residue_map.find(key) == pdb_residue_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_residue_map.find(key) == pdb_residue_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;

        models_->SetModels(models);
        this->UpdateConnectCard();
    }
}

void PdbFile::DeleteAtom(PdbFileSpace::PdbAtomCard* target_atom)
{
    std::string target_residue_name = target_atom->GetAtomResidueName();
    char target_residue_chain_id = target_atom->GetAtomChainId();
    int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
    char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
    char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_atoms[serial_number] = atom;
                        updated_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_heterogen_atoms[serial_number] = atom;
                        updated_heterogen_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
    this->UpdateConnectCard();
}

void PdbFile::DeleteAtoms(PdbAtomCardVector target_atoms)
{
    std::map<std::string, PdbFileSpace::PdbAtomCard* > pdb_atom_map = std::map<std::string, PdbFileSpace::PdbAtomCard* >();
    for(PdbAtomCardVector::iterator it = target_atoms.begin(); it != target_atoms.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* target_atom = (*it);
        std::string target_residue_name = target_atom->GetAtomResidueName();
        char target_residue_chain_id = target_atom->GetAtomChainId();
        int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
        char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
        char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        std::string target_key = ss.str();
        pdb_atom_map[target_key] = target_atom;
    }

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_atom_map.find(key) == pdb_atom_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    PdbFileSpace::PdbAtomCard* target_atom = pdb_atom_map[key];
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_atoms[serial_number] = atom;
                        updated_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_atom_map.find(key) == pdb_atom_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    PdbFileSpace::PdbAtomCard* target_atom = pdb_atom_map[key];
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_heterogen_atoms[serial_number] = atom;
                        updated_heterogen_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
    this->UpdateConnectCard();
}

void PdbFile::DeleteAtomWithTheGivenModelNumber(PdbFileSpace::PdbAtomCard* target_atom, int model_number)
{
    std::string target_residue_name = target_atom->GetAtomResidueName();
    char target_residue_chain_id = target_atom->GetAtomChainId();
    int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
    char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
    char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_atoms[serial_number] = atom;
                        updated_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_heterogen_atoms[serial_number] = atom;
                        updated_heterogen_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;
        models_->SetModels(models);
        this->UpdateConnectCard();
    }
}

void PdbFile::DeleteAtomsWithTheGivenModelNumber(PdbAtomCardVector target_atoms, int model_number)
{
    std::map<std::string, PdbFileSpace::PdbAtomCard* > pdb_atom_map = std::map<std::string, PdbFileSpace::PdbAtomCard* >();
    for(PdbAtomCardVector::iterator it = target_atoms.begin(); it != target_atoms.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* target_atom = (*it);
        std::string target_residue_name = target_atom->GetAtomResidueName();
        char target_residue_chain_id = target_atom->GetAtomChainId();
        int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
        char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
        char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        std::string target_key = ss.str();
        pdb_atom_map[target_key] = target_atom;
    }

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_atom_map.find(key) == pdb_atom_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    PdbFileSpace::PdbAtomCard* target_atom = pdb_atom_map[key];
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_atoms[serial_number] = atom;
                        updated_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_atom_map.find(key) == pdb_atom_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    PdbFileSpace::PdbAtomCard* target_atom = pdb_atom_map[key];
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_heterogen_atoms[serial_number] = atom;
                        updated_heterogen_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;
        models_->SetModels(models);
        this->UpdateConnectCard();
    }
}

void PdbFile::UpdateResidueName(PdbResidue *residue, std::string updated_residue_name)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
             char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) == 0)
                {
                    atom->SetAtomResidueName(updated_residue_name);
                    updated_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_atoms_vector.push_back(atom);
                }
                else
                {
                    updated_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_atoms_vector.push_back(atom);
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
                char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) == 0)
                {
                    atom->SetAtomResidueName(updated_residue_name);
                    updated_heterogen_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                }
                else
                {
                    updated_heterogen_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
}

void PdbFile::UpdateResidueNameWithTheGivenModelNumber(PdbResidue *residue, std::string updated_residue_name, int model_number)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
                char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) == 0)
                {
                    atom->SetAtomResidueName(updated_residue_name);
                    updated_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_atoms_vector.push_back(atom);
                }
                else
                {
                    updated_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_atoms_vector.push_back(atom);
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                //char alternate_location = atom->GetAtomAlternateLocation();
                char alternate_location = gmml::BLANK_SPACE;
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) == 0)
                {
                    atom->SetAtomResidueName(updated_residue_name);
                    updated_heterogen_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                }
                else
                {
                    updated_heterogen_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;

        models_->SetModels(models);
    }
}

void PdbFile::InsertResidueBeforeWithTheGivenModelNumber(PdbAtomSection* residue, int model_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card = new PdbAtomSection();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();
            // bool located = true; // The order is actually wrong. Can't find what orders it. I'll just look for the correct name. Hopefully we can burn this code and start again.
            bool located = false;

            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_of_residue = residue->GetOrderedAtomCards();
            PdbFileSpace::PdbAtomCard* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                if (atom->GetAtomName() == "N")
                {
                    located = true;
                }
                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                    {
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    }
                    if(located)
                    {
                        // : update coordinates with respect to it2
                        GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        GeometryTopology::Coordinate refCoordInNewResidue;
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                            if ((*it3)->GetAtomName() == "C") // Only works for ACE. This code needs to die.
                            {
                                refCoordInNewResidue  = (*it3)->GetAtomOrthogonalCoordinate();
                                gmml::log(__LINE__, __FILE__, gmml::INF, "Found the atom named C" );
                            }
                        }
                        gmml::log(__LINE__, __FILE__,  gmml::INF, "These code triggered: " + (*it2)->GetAtomName() + "_" + (*it2)->GetAtomResidueName() );
                        //base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, -1);
                        base_coordinate->TranslateCoordinateSet(coordinate_set, gmml::BOND_LENGTH, refCoordInNewResidue);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, -1);
                        //base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, -1);
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom->GetAlternateAtomCards());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                        sequence_number++;
                        offset++;
                        located = false;
                    }
                    if(!located)
                    {
                        PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(),
                                                            atom->GetAlternateAtomCards());
                        updated_atoms[serial_number] = updated_atom;
                        updated_atoms_vector.push_back(updated_atom);
                        serial_number++;
                        if(sequence_number != atom->GetAtomResidueSequenceNumber())
                            sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    }
                }
                else
                {
                    if(it2 != ordered_atoms.begin())
                    {

                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    updated_atoms[serial_number] = updated_atom;
                    updated_atoms_vector.push_back(updated_atom);
                    serial_number++;
                    if(sequence_number != atom->GetAtomResidueSequenceNumber())
                        sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                }
            }
            updated_atom_card->SetAtomCards(updated_atoms);
            updated_atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection* updated_heterogen_atom_card = new PdbHeterogenAtomSection();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector  updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2);
                PdbFileSpace::PdbAtomCard* updated_heterogen_atom = new PdbFileSpace::PdbAtomCard(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge(),
                                                              heterogen_atom->GetAlternateAtomCards());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
        models_->SetModels(updated_models);
    }
}

// OG Dec2021: This code is nuts. It's getting ordered atoms for some reason, but they aren't ordered correctly so then translation happens to the wrong atom. There's this located bool that's mental.
// I'm adding in a shortcut that just looks for the name, this whole code/approach needs to go away.
void PdbFile::InsertResidueAfterWithTheGivenModelNumber(PdbAtomSection* residue, int model_number)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "InsertResidueAfterWithTheGivenModelNumber");
    std::map<int, PdbModelCard*> modelNumberToModelCardMap = models_->GetModels();
    if(modelNumberToModelCardMap.empty())
    {
        return; // Better to throw, just mirroring actual behaviour for now.
    }
    PdbModelCard* model = modelNumberToModelCardMap[model_number];
    PdbModelResidueSet::AtomCardVector atom_cards = model->GetModelResidueSet()->GetAtomCards();
    PdbModelResidueSet::AtomCardVector updated_atom_cards;
    int serial_number = 1;
    int sequence_number = 1;
    int offset = 0;
    PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_of_residue = residue->GetOrderedAtomCards();
     PdbFileSpace::PdbAtomCard* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
     char residue_chain_id = first_atom_in_residue->GetAtomChainId();
     int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection* updated_atom_card = new PdbAtomSection();
        updated_atom_card->SetRecordName(atom_card->GetRecordName());
        PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();

        PdbAtomSection::PdbAtomMap updated_atoms;
        PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();
        bool located = false;
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
            if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "residue_chain_id is " + residue_chain_id);
                gmml::log(__LINE__, __FILE__, gmml::INF, "residue_sequence_number " + residue_sequence_number);
                if(it2 != ordered_atoms.begin())
                {
                    PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                    it2++;
                    if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                    {
                        int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                        if(diff == 1)
                            sequence_number++;
                        else
                            sequence_number += diff;
                    }
                }
                else
                {
                    sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                }
                PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                updated_atoms[serial_number] = updated_atom;
                updated_atoms_vector.push_back(updated_atom);
                serial_number++;
                if(sequence_number != atom->GetAtomResidueSequenceNumber())
                    sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                located = true;
            }
            else
            {
                if(located)
                {
                    //: update coordinates with respect to it2
//                    GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
//                    GeometryTopology::Coordinate base_coordinate = proteinCAtom.GetAtomOrthogonalCoordinate();
//                    GeometryTopology::Coordinate coordOfNAtomInNME, coordOfCH3AtomInNME, coordOfHAtomInNME;
//                    for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
//                    {
//                        coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
//                        if ((*it3)->GetAtomName() == "N") // Only works for NME. This code needs to die.
//                        {
//                            coordOfNAtomInNME  = (*it3)->GetAtomOrthogonalCoordinate();
//                            gmml::log(__LINE__, __FILE__, gmml::INF, "Found the atom named N in NME" );
//                        }
//                        if ((*it3)->GetAtomName() == "CH3") // Only works for NME. This code needs to die.
//                        {
//                            coordOfCH3AtomInNME  = (*it3)->GetAtomOrthogonalCoordinate();
//                            gmml::log(__LINE__, __FILE__, gmml::INF, "Found the atom named CH3 in NME" );
//                        }
//                        if ((*it3)->GetAtomName() == "H") // Only works for NME. This code needs to die.
//                        {
//                            coordOfHAtomInNME  = (*it3)->GetAtomOrthogonalCoordinate();
//                            gmml::log(__LINE__, __FILE__, gmml::INF, "Found the atom named H in NME" );
//                        }
//                    }
//                    gmml::log(__LINE__, __FILE__,  gmml::INF, "Thon code triggered: " + proteinCAtom.GetAtomName() + "_" + proteinCAtom.GetAtomResidueName() );
                    //base_coordinate->TranslateCoordinateSet(coordinate_set, gmml::BOND_LENGTH, coordOfNAtomInNME);
                   // base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                  //  base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
//                    base_coordinate.TranslateCoordinateSet(coordinate_set, gmml::BOND_LENGTH, coordOfNAtomInNME);
//                    base_coordinate.SetAngle(coordinate_set, 180.0, coordOfNAtomInNME, coordOfCH3AtomInNME);
//                    GeometryTopology::Coordinate proteinOAtomCoord = proteinOAtom.GetAtomOrthogonalCoordinate();
//                    GeometryTopology::SetDihedralAngle(180.0, &proteinOAtomCoord, &base_coordinate, coordOfNAtomInNME, coordOfHAtomInNME, coordinate_set);
////                    proteinOAtom.GetAtomOrthogonalCoordinate().SetDihedralAngle(coordinate_set, 180.0, base_coordinate, coordOfNAtomInNME, coordOfHAtomInNME);
                    sequence_number++;
                    for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                    {
                        PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                        int index = distance(ordered_atoms_of_residue.begin(), it3);
                        PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom->GetAlternateAtomCards());
                        updated_atoms[serial_number] = new_atom;
                        updated_atoms_vector.push_back(new_atom);
                        serial_number++;
                    }
                    sequence_number++;
                    offset++;
                    located = false;
                }
                if(!located)
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                    {
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    }
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(),
                            atom->GetAlternateAtomCards());
                    updated_atoms[serial_number] = updated_atom;
                    updated_atoms_vector.push_back(updated_atom);
                    serial_number++;
                    if(sequence_number != atom->GetAtomResidueSequenceNumber())
                        sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                }
            }
            if(it2 == --ordered_atoms.end() && located)
            {
                // : update coordinates with respect to it2
                GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                GeometryTopology::CoordinateVector justATestCoordinateSetBro;
                GeometryTopology::Coordinate base_coordinate = proteinCAtom.GetAtomOrthogonalCoordinate();
                GeometryTopology::Coordinate *coordOfNAtomInNME, *coordOfCH3AtomInNME, *coordOfHAtomInNME;
                for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                {
                    coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                    if ((*it3)->GetAtomName() == "N") // Only works for NME. This code needs to die.
                    {
                        coordOfNAtomInNME = new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate());
                        gmml::log(__LINE__, __FILE__, gmml::INF, "Found the atom named N in NME" );
                    }
                    else if ((*it3)->GetAtomName() == "CH3") // Only works for NME. This code needs to die.
                    {
                        coordOfCH3AtomInNME = new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate());
                        gmml::log(__LINE__, __FILE__, gmml::INF, "Found the atom named CH3 in NME" );
                    }
                    else if ((*it3)->GetAtomName() == "H") // Only works for NME. This code needs to die.
                    {
                        coordOfHAtomInNME = new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate());
                        gmml::log(__LINE__, __FILE__, gmml::INF, "Found the atom named H in NME" );
                    }
                    else
                    {
                        justATestCoordinateSetBro.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                    }
                }
                base_coordinate.TranslateCoordinateSet(coordinate_set, gmml::BOND_LENGTH, coordOfNAtomInNME);
                //base_coordinate.SetAngle(coordinate_set, 125.0, coordOfNAtomInNME, coordOfCH3AtomInNME);
                //proteinOAtom.GetAtomOrthogonalCoordinate().SetDihedralAngle(justATestCoordinateSetBro, 180.0, base_coordinate, coordOfNAtomInNME, coordOfHAtomInNME);
                GeometryTopology::Coordinate proteinOAtomCoord = proteinOAtom.GetAtomOrthogonalCoordinate();

                //GeometryTopology::SetDihedralAngle(180.0, &proteinOAtomCoord, &base_coordinate, coordOfNAtomInNME, coordOfHAtomInNME, coordinate_set);

//                base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                //base_coordinate.RotateAngularAll(coordinate_set, 125.0, 1);
            //    base_coordinate.RotateAngularAll(CoordinateVector coordinate_set, 125, refCoordInNewResidue, Coordinate* a3);

                gmml::log(__LINE__, __FILE__,  gmml::INF, "That code triggered: " + proteinCAtom.GetAtomName() + "_" + proteinCAtom.GetAtomResidueName() );
                //base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                sequence_number++;
                for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                {
                    PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                    int index = distance(ordered_atoms_of_residue.begin(), it3);
                    PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom->GetAlternateAtomCards());
                    updated_atoms[serial_number] = new_atom;
                    updated_atoms_vector.push_back(new_atom);
                    serial_number++;
                }
                sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                located = false;
                sequence_number++;
                offset++;
            }
        }
        updated_atom_card->SetAtomCards(updated_atoms);
        updated_atom_card->SetOrderedAtomCards(updated_atoms_vector);
        updated_atom_cards.push_back(updated_atom_card);
        serial_number++;
    }
    PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
    for(auto &heterogen_atom_card : model->GetModelResidueSet()->GetHeterogenAtomCards())
    {
        PdbHeterogenAtomSection* updated_heterogen_atom_card = new PdbHeterogenAtomSection();
        updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
        PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
        PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
        PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
        for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2);
            PdbFileSpace::PdbAtomCard* updated_heterogen_atom = new PdbFileSpace::PdbAtomCard(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                    heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                    heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                    heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                    heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge(),
                    heterogen_atom->GetAlternateAtomCards());
            updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
            updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
            serial_number++;
        }
        updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
        updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
        updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
    }
    PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
    updated_residue_set->SetAtomCards(updated_atom_cards);
    updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
    PdbModelCard* updated_model = new PdbModelCard();
    updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
    updated_model->SetModelResidueSet(updated_residue_set);
    PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
    updated_models[updated_model->GetModelSerialNumber()] = updated_model;
    models_->SetModels(updated_models);
}

void PdbFile::SplitAtomCardOfModelCard(char split_point_chain_id, int split_point_sequence_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card_first_part = new PdbAtomSection();
            PdbAtomSection* updated_atom_card_second_part = new PdbAtomSection();
            updated_atom_card_first_part->SetRecordName(atom_card->GetRecordName());
            updated_atom_card_second_part->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap atoms_first_part;
            PdbAtomSection::PdbAtomMap atoms_second_part;
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_first_part = PdbAtomSection::PdbAtomCardOrderVector();
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_second_part = PdbAtomSection::PdbAtomCardOrderVector();

            bool located = false;
            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                char residue_chain_id = atom->GetAtomChainId();
                int residue_sequence_number = atom->GetAtomResidueSequenceNumber();

                if(residue_chain_id == split_point_chain_id && residue_sequence_number == split_point_sequence_number)
                {
                    located = true;
                }
                if(located)
                {
                    serial_number++;
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    atoms_second_part[serial_number] = updated_atom;
                    ordered_atoms_second_part.push_back(updated_atom);
                }
                if(!located)
                {
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    atoms_first_part[serial_number] = updated_atom;
                    ordered_atoms_first_part.push_back(updated_atom);
                    serial_number++;
                }
            }
            updated_atom_card_first_part->SetAtomCards(atoms_first_part);
            updated_atom_card_first_part->SetOrderedAtomCards(ordered_atoms_first_part);
            updated_atom_card_second_part->SetAtomCards(atoms_second_part);
            updated_atom_card_second_part->SetOrderedAtomCards(ordered_atoms_second_part);
            if(ordered_atoms_first_part.size() > 0)
            {
                updated_atom_cards.push_back(updated_atom_card_first_part);
                serial_number++;
            }
            if(ordered_atoms_second_part.size() > 0)
            {
                updated_atom_cards.push_back(updated_atom_card_second_part);
                serial_number++;
            }
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();

        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
    }
    models_->SetModels(updated_models);
}

void PdbFile::SplitAtomCardOfModelCardWithTheGivenModelNumber(char split_point_chain_id, int split_point_sequence_number, int model_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
        PdbModelCard* model = models[model_number];
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card_first_part = new PdbAtomSection();
            PdbAtomSection* updated_atom_card_second_part = new PdbAtomSection();
            updated_atom_card_first_part->SetRecordName(atom_card->GetRecordName());
            updated_atom_card_second_part->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap atoms_first_part;
            PdbAtomSection::PdbAtomMap atoms_second_part;
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_first_part = PdbAtomSection::PdbAtomCardOrderVector();
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_second_part = PdbAtomSection::PdbAtomCardOrderVector();

            bool located = false;
            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                char residue_chain_id = atom->GetAtomChainId();
                int residue_sequence_number = atom->GetAtomResidueSequenceNumber();

                if(residue_chain_id == split_point_chain_id && residue_sequence_number == split_point_sequence_number)
                {
                    located = true;
                }
                if(located)
                {
                    serial_number++;
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    atoms_second_part[serial_number] = updated_atom;
                    ordered_atoms_second_part.push_back(updated_atom);
                }
                if(!located)
                {
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    atoms_first_part[serial_number] = updated_atom;
                    ordered_atoms_first_part.push_back(updated_atom);
                    serial_number++;
                }
            }
            updated_atom_card_first_part->SetAtomCards(atoms_first_part);
            updated_atom_card_first_part->SetOrderedAtomCards(ordered_atoms_first_part);
            updated_atom_card_second_part->SetAtomCards(atoms_second_part);
            updated_atom_card_second_part->SetOrderedAtomCards(ordered_atoms_second_part);
            if(ordered_atoms_first_part.size() > 0)
            {
                updated_atom_cards.push_back(updated_atom_card_first_part);
                serial_number++;
            }
            if(ordered_atoms_second_part.size() > 0)
            {
                updated_atom_cards.push_back(updated_atom_card_second_part);
                serial_number++;
            }
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();

        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
        models_->SetModels(updated_models);
    }
}
