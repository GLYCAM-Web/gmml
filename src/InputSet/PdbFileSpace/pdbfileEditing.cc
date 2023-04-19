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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards = PdbModelResidueSet::PdbAtomSectionVector();
        int serial_number = 1;
        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomSectionVector();
        for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
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
    gmml::log(__LINE__, __FILE__, gmml::INF, "InsertResidueBeforeWithTheGivenModelNumber() aka ACEInserter was called");
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        PdbFileSpace::PdbAtomCard *cAtomInProtein, *nAtomInProtein, *caAtomInProtein; // OG Dec 21 Sorry about the warnings.
        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card = new PdbAtomSection();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();
            bool located = true;

            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_of_residue = residue->GetOrderedAtomCards();
            PdbFileSpace::PdbAtomCard* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    // OG Go forward from here and get all the info you need.
                    PdbAtomSection::PdbAtomCardOrderVector::iterator fromHere = it2;
                    while((fromHere != ordered_atoms.end()) && ((*fromHere)->GetAtomResidueSequenceNumber() == residue_sequence_number) && (located == true))
                    {
                        PdbFileSpace::PdbAtomCard* futureAtom = (*fromHere);
                        std::stringstream ss;
                        ss << residue_chain_id << "_" << residue_sequence_number << "_" << atom->GetAtomResidueName();
                        if (futureAtom->GetAtomName() == "C") // This code needs to die.
                        {
                            cAtomInProtein = futureAtom;
                            gmml::log(__LINE__, __FILE__, gmml::INF, "ACE. Found the atom named C in " + ss.str());
                        }
                        if (futureAtom->GetAtomName() == "N") // This code needs to die.
                        {
                            nAtomInProtein = futureAtom;
                            gmml::log(__LINE__, __FILE__, gmml::INF, "ACE. Found the atom named N in " + ss.str());
                        }
                        if (futureAtom->GetAtomName() == "CA") // This code needs to die.
                        {
                            caAtomInProtein = futureAtom;
                            gmml::log(__LINE__, __FILE__, gmml::INF, "ACE. Found the atom named CA in " + ss.str());
                        }
                        fromHere++;
                    }
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
                    if(located)
                    {
                        sequence_number = atom->GetAtomResidueSequenceNumber() - 1; // OG edit: Single gaps will end up with the same ACE NME resid numbers. Otherwise good.
                        GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                        GeometryTopology::Coordinate* cCoordProtein = new GeometryTopology::Coordinate(cAtomInProtein->GetAtomOrthogonalCoordinate());
                        GeometryTopology::Coordinate* caCoordProtein = new GeometryTopology::Coordinate(caAtomInProtein->GetAtomOrthogonalCoordinate());
                        GeometryTopology::Coordinate* nCoordProtein = new GeometryTopology::Coordinate(nAtomInProtein->GetAtomOrthogonalCoordinate());
                        // This is a bad idea, but it works. Better to load in the template and align it, but the support code didn't work and I'm a quitter.
                        GeometryTopology::Coordinate cCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(cCoordProtein, caCoordProtein, nCoordProtein, 120.0, -130.0, 1.4);
                        GeometryTopology::Coordinate oCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, cCoordACE, 120.0, 0.0, 1.23);
                        GeometryTopology::Coordinate ch3CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, cCoordACE, 125.0, 180.0, 1.48);
                        GeometryTopology::Coordinate hh31CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 180.0, 1.09);
                        GeometryTopology::Coordinate hh32CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 60.0, 1.09);
                        GeometryTopology::Coordinate hh33CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, -60.0, 1.09);
                        // Quantity and order matter.
                        coordinate_set.push_back(new GeometryTopology::Coordinate(hh31CoordACE));
                        coordinate_set.push_back(new GeometryTopology::Coordinate(ch3CoordACE));
                        coordinate_set.push_back(new GeometryTopology::Coordinate(hh32CoordACE));
                        coordinate_set.push_back(new GeometryTopology::Coordinate(hh33CoordACE));
                        coordinate_set.push_back(new GeometryTopology::Coordinate(cCoordACE));
                        coordinate_set.push_back(new GeometryTopology::Coordinate(oCoordACE));
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number,
                                                            atom_of_residue->GetAtomName(),
                                                            atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),
                                                            atom_of_residue->GetAtomChainId(),
                                                            sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(),
                                                            coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(),
                                                            atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(),
                                                            atom_of_residue->GetAtomCharge(),
                                                            atom->GetAlternateAtomCards());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        gmml::log(__LINE__, __FILE__, gmml::INF, "Created ACE residue numbered: " + std::to_string(sequence_number));
                        sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                        sequence_number++;
                       // offset++;
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
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
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;

        models_->SetModels(updated_models);
    }
}

void PdbFile::InsertResidueAfterWithTheGivenModelNumber(PdbAtomSection* residue, int model_number)
{
    std::stringstream ss;
    ss << residue->GetRecordName();
    gmml::log(__LINE__, __FILE__, gmml::INF, "InsertResidueAfterWithTheGivenModelNumber() aka NME inserter was called.");
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() == 0)
    {
        return;
    }
    PdbModelCard* model = models[model_number];
    PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
    PdbModelCard* updated_model = new PdbModelCard();
    updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
    PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
    PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
    int serial_number = 1;
    int sequence_number = 1;
    int offset = 0;
    PdbFileSpace::PdbAtomCard *cAtomInProtein, *oAtomInProtein, *caAtomInProtein; // OG Dec 21 Sorry about the warnings.
    for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection* updated_atom_card = new PdbAtomSection();
        updated_atom_card->SetRecordName(atom_card->GetRecordName());
        PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
        PdbAtomSection::PdbAtomMap updated_atoms;
        PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();
        bool located = false;

        PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_of_residue = residue->GetOrderedAtomCards();
        PdbFileSpace::PdbAtomCard* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
        char residue_chain_id = first_atom_in_residue->GetAtomChainId();
        int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();
        std::string residue_name = first_atom_in_residue->GetAtomResidueName();
        std::stringstream ss1;
        ss1 << residue_chain_id << "_" << residue_sequence_number << "_" << residue_name;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Residue passed in is: " + ss1.str() );
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);

            if((residue_chain_id == atom->GetAtomChainId()) && (residue_sequence_number == atom->GetAtomResidueSequenceNumber()) && (atom->GetAtomResidueName() != "ACE") ) // this code is making me commit crimes against humanity.
            {
                std::stringstream ss;
                ss << residue_chain_id << "_" << residue_sequence_number << "_" << atom->GetAtomResidueName();
                if (atom->GetAtomName() == "C") // This code needs to die.
                {
                    cAtomInProtein = atom;
                    gmml::log(__LINE__, __FILE__, gmml::INF, "NME: Found the atom named C in " + ss.str());
                }
                if (atom->GetAtomName() == "O") // This code needs to die.
                {
                    oAtomInProtein = atom;
                    gmml::log(__LINE__, __FILE__, gmml::INF, "NME: Found the atom named O in " + ss.str());
                }
                if (atom->GetAtomName() == "CA") // This code needs to die.
                {
                    caAtomInProtein = atom;
                    gmml::log(__LINE__, __FILE__, gmml::INF, "NME: Found the atom named CA in " + ss.str());
                }
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
                    GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector(); // These are pointers to newly created copies!
                    GeometryTopology::Coordinate* cAtomInProteinCoord = new GeometryTopology::Coordinate(cAtomInProtein->GetAtomOrthogonalCoordinate());
                    GeometryTopology::Coordinate* caAtomInProteinCoord = new GeometryTopology::Coordinate(caAtomInProtein->GetAtomOrthogonalCoordinate());
                    GeometryTopology::Coordinate* oAtomInProteinCoord = new GeometryTopology::Coordinate(oAtomInProtein->GetAtomOrthogonalCoordinate());
                    // This is a bad idea, but it works. Better to load in the template and align it, but the support code didn't work.
                    GeometryTopology::Coordinate bcoordOfNAtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(oAtomInProteinCoord, caAtomInProteinCoord, cAtomInProteinCoord, 120.0, 180.0, 1.4);
                    coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfNAtomInNME));
                    GeometryTopology::Coordinate bcoordOfHAtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(oAtomInProteinCoord, cAtomInProteinCoord, bcoordOfNAtomInNME, 109.0, 180.0, 1.0);
                    coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfHAtomInNME));
                    GeometryTopology::Coordinate bcoordOfCH3AtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(caAtomInProteinCoord, cAtomInProteinCoord, bcoordOfNAtomInNME, 125.0, 180.0, 1.48);
                    coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfCH3AtomInNME));
                    GeometryTopology::Coordinate bcoordOfHH31AtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(bcoordOfHAtomInNME, bcoordOfNAtomInNME, bcoordOfCH3AtomInNME, 109.0, 180.0, 1.09);
                    coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfHH31AtomInNME));
                    GeometryTopology::Coordinate bcoordOfHH32AtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(bcoordOfHAtomInNME, bcoordOfNAtomInNME, bcoordOfCH3AtomInNME, 109.0, 60.0, 1.09);
                    coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfHH32AtomInNME));
                    GeometryTopology::Coordinate bcoordOfHH33AtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(bcoordOfHAtomInNME, bcoordOfNAtomInNME, bcoordOfCH3AtomInNME, 109.0, -60.0, 1.09);
                    coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfHH33AtomInNME));
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
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Created NME residue numbered: " + std::to_string(sequence_number));
                    sequence_number++;
                   // offset++;
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
            if(it2 == --ordered_atoms.end() && located) // Note that the below no longer uses it2, I've made the spaghetti worse.
            {
                GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector(); // These are pointers to newly created copies!
                GeometryTopology::Coordinate* cAtomInProteinCoord = new GeometryTopology::Coordinate(cAtomInProtein->GetAtomOrthogonalCoordinate());
                GeometryTopology::Coordinate* caAtomInProteinCoord = new GeometryTopology::Coordinate(caAtomInProtein->GetAtomOrthogonalCoordinate());
                GeometryTopology::Coordinate* nAtomInProteinCoord = new GeometryTopology::Coordinate(oAtomInProtein->GetAtomOrthogonalCoordinate());
                // This is a bad idea, but it works. Better to load in the template and align it, but the support code didn't work.
                GeometryTopology::Coordinate bcoordOfNAtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(nAtomInProteinCoord, caAtomInProteinCoord, cAtomInProteinCoord, 120.0, 140.0, 1.4);
                coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfNAtomInNME));
                GeometryTopology::Coordinate bcoordOfHAtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(caAtomInProteinCoord, cAtomInProteinCoord, bcoordOfNAtomInNME, 109.0, 0.0, 1.0);
                coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfHAtomInNME));
                GeometryTopology::Coordinate bcoordOfCH3AtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(caAtomInProteinCoord, cAtomInProteinCoord, bcoordOfNAtomInNME, 125.0, 180.0, 1.48);
                coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfCH3AtomInNME));
                GeometryTopology::Coordinate bcoordOfHH31AtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(bcoordOfHAtomInNME, bcoordOfNAtomInNME, bcoordOfCH3AtomInNME, 109.0, 180.0, 1.09);
                coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfHH31AtomInNME));
                GeometryTopology::Coordinate bcoordOfHH32AtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(bcoordOfHAtomInNME, bcoordOfNAtomInNME, bcoordOfCH3AtomInNME, 109.0, 60.0, 1.09);
                coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfHH32AtomInNME));
                GeometryTopology::Coordinate bcoordOfHH33AtomInNME = GeometryTopology::get_cartesian_point_from_internal_coords(bcoordOfHAtomInNME, bcoordOfNAtomInNME, bcoordOfCH3AtomInNME, 109.0, -60.0, 1.09);
                coordinate_set.push_back(new GeometryTopology::Coordinate(bcoordOfHH33AtomInNME));
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
                gmml::log(__LINE__, __FILE__, gmml::INF, "Created NME residue numbered: " + std::to_string(sequence_number));
                sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                located = false;
                sequence_number++;
                //offset++;
            }
        }
        updated_atom_card->SetAtomCards(updated_atoms);
        updated_atom_card->SetOrderedAtomCards(updated_atoms_vector);
        updated_atom_cards.push_back(updated_atom_card);
        serial_number++;
    }
    PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
    PdbModelResidueSet::HeterogenAtomSectionVector updated_heterogen_atom_cards;
    for(PdbModelResidueSet::HeterogenAtomSectionVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
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
    updated_residue_set->SetAtomCards(updated_atom_cards);
    updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
    updated_model->SetModelResidueSet(updated_residue_set);
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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        int serial_number = 1;
        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();

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
        PdbModelResidueSet::PdbAtomSectionVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::PdbAtomSectionVector updated_atom_cards;
        int serial_number = 1;
        for(PdbModelResidueSet::PdbAtomSectionVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
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
        PdbModelResidueSet::HeterogenAtomSectionVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();

        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
        models_->SetModels(updated_models);
    }
}
