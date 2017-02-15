// Author: Alireza Khatamian

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <exception>
#include <cctype>

#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheadercard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcompoundcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodeltypecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequence.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodification.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogencard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogen.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenname.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonym.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbformulacard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbformula.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbhelixcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbhelix.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheet.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfidebondcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblink.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsitecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsite.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcrystallographiccard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdboriginxncard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdboriginxn.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbscalencard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbscalen.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmatrixn.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresidue.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/coordinate.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFile::PdbFile()
{
    path_ = "GMML-Generated";
    header_ = NULL;
    title_ = NULL;
    compound_ = NULL;
    number_of_models_ = NULL;
    model_type_ = NULL;
    residues_sequence_ = NULL;
    residue_modification_ = NULL;
    heterogens_ = NULL;
    heterogens_name_ = NULL;
    heterogen_synonyms_ = NULL;
    formulas_ = NULL;
    helixes_ = NULL;
    sheets_ = NULL;
    disulfide_bonds_ = NULL;
    links_ = NULL;
    sites_ = NULL;
    crystallography_ = NULL;
    origins_ = NULL;
    scales_ = NULL;
    matrices_ = NULL;
    models_ = NULL;
    connectivities_ = NULL;
    serial_number_mapping_ = PdbFile::PdbSerialNumberMapping();
    sequence_number_mapping_ = PdbFile::PdbSequenceNumberMapping();
}

PdbFile::PdbFile(const std::string &pdb_file)
{
    path_ = pdb_file;
    header_ = NULL;
    title_ = NULL;
    compound_ = NULL;
    number_of_models_ = NULL;
    model_type_ = NULL;
    residues_sequence_ = NULL;
    residue_modification_ = NULL;
    heterogens_ = NULL;
    heterogens_name_ = NULL;
    heterogen_synonyms_ = NULL;
    formulas_ = NULL;
    helixes_ = NULL;
    sheets_ = NULL;
    disulfide_bonds_ = NULL;
    links_ = NULL;
    sites_ = NULL;
    crystallography_ = NULL;
    origins_ = NULL;
    scales_ = NULL;
    matrices_ = NULL;
    models_ = NULL;
    connectivities_ = NULL;
    serial_number_mapping_ = PdbFile::PdbSerialNumberMapping();
    sequence_number_mapping_ = PdbFile::PdbSequenceNumberMapping();

    std::ifstream in_file;
    if(std::ifstream(pdb_file.c_str()))
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Opening PDB file ...");
        cout << "Opening PDB file ..." << endl;
        in_file.open(pdb_file.c_str());
    }
    else
    {
        throw PdbFileProcessingException(__LINE__, "PDB file not found");
    }

    string line = "";
    string temp = "";
    stringstream ss;
    while(!in_file.eof())
    {
        if(!getline(in_file, line))
            break;
        else
        {
            temp = line.substr(0,6);
            temp = Trim(temp);
            if(temp.find("END") != string::npos || temp.compare("END") == 0)
                break;
            else if(!line.empty())
                ss << line << endl;
        }
    }
    in_file.close();
    if(temp.find("END") == string::npos || temp.compare("END") != 0)
    {
        std::ofstream out_file;
        out_file.open(pdb_file.c_str());
        out_file << ss.str() << "END";
        out_file.close();
    }
    else
    {
        std::ofstream out_file;
        out_file.open(pdb_file.c_str());
        out_file << ss.str() << temp;
        out_file.close();
    }
    in_file.open(pdb_file.c_str());
    if(!Read(in_file))
    {
        throw PdbFileProcessingException(__LINE__, "Reading PDB file exception");
    }
    in_file.close();            /// Close the pdb files
}
PdbFile* PdbFile::LoadPdbFile()
{
    PdbFile* pdb = new PdbFile();
    return pdb;
}
PdbFile* PdbFile::LoadPdbFile(const std::string &pdb_file)
{
    PdbFile* pdb = new PdbFile(pdb_file);
    return pdb;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbFile::GetPath()
{
    return path_;
}

PdbHeaderCard* PdbFile::GetHeader()
{
    return header_;
}

PdbTitleCard* PdbFile::GetTitle()
{
    return title_;
}

PdbCompoundCard* PdbFile::GetCompound()
{
    return compound_;
}

PdbNumModelCard* PdbFile::GetNumberOfModels()
{
    return number_of_models_;
}

PdbModelTypeCard* PdbFile::GetModelType()
{
    return model_type_;
}

PdbResidueSequenceCard* PdbFile::GetResiduesSequence()
{
    return residues_sequence_;
}

PdbResidueModificationCard* PdbFile::GetResidueModification()
{
    return residue_modification_;
}

PdbHeterogenCard* PdbFile::GetHeterogens()
{
    return heterogens_;
}

PdbHeterogenNameCard* PdbFile::GetHeterogensName()
{
    return heterogens_name_;
}

PdbHeterogenSynonymCard* PdbFile::GetHeterogenSynonyms()
{
    return heterogen_synonyms_;
}

PdbFormulaCard* PdbFile::GetFormulas()
{
    return formulas_;
}

PdbHelixCard* PdbFile::GetHelixes()
{
    return helixes_;
}

PdbSheetCard* PdbFile::GetSheets()
{
    return sheets_;
}

PdbDisulfideBondCard* PdbFile::GetDisulfideBonds()
{
    return disulfide_bonds_;
}

PdbLinkCard* PdbFile::GetLinks()
{
    return links_;
}

PdbSiteCard* PdbFile::GetSites()
{
    return sites_;
}

PdbCrystallographicCard* PdbFile::GetCrystallography()
{
    return crystallography_;
}

PdbOriginXnCard* PdbFile::GetOrigins()
{
    return origins_;
}

PdbScaleNCard* PdbFile::GetScales()
{
    return scales_;
}

PdbMatrixNCard* PdbFile::GetMatrices()
{
    return matrices_;
}

PdbModelCard* PdbFile::GetModels()
{
    return models_;
}

PdbConnectCard* PdbFile::GetConnectivities()
{
    return connectivities_;
}
PdbFile::PdbSerialNumberMapping PdbFile::GetSerialNumberMapping()
{
    return serial_number_mapping_;
}
PdbFile::PdbSerialNumberMapping PdbFile::GetSequenceNumberMapping()
{
    return sequence_number_mapping_;
}

PdbFile::PdbPairVectorAtomNamePositionFlag PdbFile::GetAllResidueNames()
{    
    PdbPairVectorAtomNamePositionFlag residue_names;
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomCard* atom_card = (*it1);
        PdbAtomCard::PdbAtomOrderVector atoms = atom_card->GetOrderedAtoms();
        for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            int dist = distance(atoms.begin(), it2);
            string atom_residue_name = atom->GetAtomResidueName();
            if(dist == 0)
            {
                pair<string, string> pair_residue_position = make_pair(atom_residue_name, "S");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
            else if(dist == atoms.size() - 1)
            {
                pair<string, string> pair_residue_position = make_pair(atom_residue_name, "E");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
            pair<string, string> pair_residue_position = make_pair(atom_residue_name, " ");
            if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                residue_names.push_back(pair_residue_position);
        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
        PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
        for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string atom_residue_name = atom->GetAtomResidueName();
            pair<string, string> pair_residue_position = make_pair(atom_residue_name, " ");
            if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
            {
                residue_names.push_back(pair_residue_position);
            }
        }
    }
    return residue_names;
}

PdbFile::PdbPairVectorAtomNamePositionFlag PdbFile::GetAllResidueNamesFromAtomCard()
{
    PdbPairVectorAtomNamePositionFlag residue_names;
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomCard* atom_card = (*it1);
        PdbAtomCard::PdbAtomOrderVector atoms = atom_card->GetOrderedAtoms();
        for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            int dist = distance(atoms.begin(), it2);
            string atom_residue_name = atom->GetAtomResidueName();
            if(dist == 0)
            {
                pair<string, string> pair_residue_position = make_pair(atom_residue_name, "S");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
            else if(dist == atoms.size() - 1)
            {
                pair<string, string> pair_residue_position = make_pair(atom_residue_name, "E");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
            else
            {
                pair<string, string> pair_residue_position = make_pair(atom_residue_name, " ");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
        }
    }
    return residue_names;
    
}

PdbFile::PdbResidueVector PdbFile::GetAllResidues()
{
    PdbFile::PdbResidueVector residues;
    map<string, bool> inserted_residues;
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomCard* atom_card = (*it1);
        PdbAtomCard::PdbAtomOrderVector atoms = atom_card->GetOrderedAtoms();
        for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            if(!inserted_residues[key])
            {
                PdbResidue* res = new PdbResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_location);
                residues.push_back(res);
                inserted_residues[key] = true;
            }
        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
        PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
        for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            if(!inserted_residues[key])
            {
                PdbResidue* res = new PdbResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_location);
                residues.push_back(res);
                inserted_residues[key] = true;
            }
        }
    }
    return residues;
}

PdbFile::PdbResidueVector PdbFile::GetAllResiduesFromAtomCard()
{
    PdbFile::PdbResidueVector residues = PdbFile::PdbResidueVector();
    map<string, bool> inserted_residues = map<string, bool>();
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomCard* atom_card = (*it1);
        PdbAtomCard::PdbAtomOrderVector atoms = atom_card->GetOrderedAtoms();
        for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            if(!inserted_residues[key])
            {
                PdbResidue* res = new PdbResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_location);
                residues.push_back(res);
                inserted_residues[key] = true;
            }
        }
    }
    return residues;
}

PdbFile::PdbAtomVector PdbFile::GetAllAtomsOfResidue(PdbResidue *residue)
{
    string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();
    
    PdbAtomVector atoms_of_residue;
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomCard* atom_card = (*it1);
        PdbAtomCard::PdbAtomOrderVector atoms = atom_card->GetOrderedAtoms();
        for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream sss;
            sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = sss.str();
            if(target_key.compare(key) == 0)
            {
                atoms_of_residue.push_back(atom);
            }
        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
        PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
        for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream sss;
            sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = sss.str();
            if(target_key.compare(key) == 0)
            {
                atoms_of_residue.push_back(atom);
            }
        }
    }
    return atoms_of_residue;
}

PdbFile::PdbResidueAtomsMap PdbFile::GetAllAtomsOfResidues()
{
    PdbFile::PdbResidueAtomsMap residue_atom_map;
    map<string, bool> inserted_residues;
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomCard* atom_card = (*it1);
        PdbAtomCard::PdbAtomOrderVector atoms = atom_card->GetOrderedAtoms();
        for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            if(!inserted_residues[key])
            {
                residue_atom_map[key] = new vector<PdbAtom*>();
                inserted_residues[key] = true;
            }
            residue_atom_map[key]->push_back(atom);
            
        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
        PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
        for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            if(!inserted_residues[key])
            {
                residue_atom_map[key] = new vector<PdbAtom*>();
                inserted_residues[key] = true;
            }
            residue_atom_map[key]->push_back(atom);
        }
    }
    return residue_atom_map;
}

PdbFile::PdbResidueAtomsMap PdbFile::GetAllAtomsInOrder(vector<string>& key_order)
{
    PdbFile::PdbResidueAtomsMap residue_atom_map;
    map<string, bool> inserted_residues;
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomCard* atom_card = (*it1);
        PdbAtomCard::PdbAtomOrderVector atoms = atom_card->GetOrderedAtoms();
        for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            if(!inserted_residues[key])
            {
                residue_atom_map[key] = new vector<PdbAtom*>();
                inserted_residues[key] = true;
                key_order.push_back(key);
            }
            residue_atom_map[key]->push_back(atom);

        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
        PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
        for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbAtom* atom = (*it2);
            string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            if(!inserted_residues[key])
            {
                residue_atom_map[key] = new vector<PdbAtom*>();
                inserted_residues[key] = true;
                key_order.push_back(key);
            }
            residue_atom_map[key]->push_back(atom);
        }
    }
    return residue_atom_map;
}

PdbFileSpace::PdbAtom* PdbFile::GetAtomOfResidueByName(PdbResidue *residue, string atom_name, PdbFile::PdbResidueAtomsMap residue_atom_map)
{
    string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();
    PdbAtomVector* atoms = residue_atom_map[target_key];
    
    for(PdbAtomVector::iterator it = atoms->begin(); it != atoms->end(); it++)
    {
        PdbAtom* atom = (*it);
        if(atom->GetAtomName().compare(atom_name) == 0)
            return atom;
    }
    return NULL;
}

PdbFileSpace::PdbAtom* PdbFile::GetAtomOfResidueByName(PdbResidue *residue, string atom_name)
{
    PdbAtomVector atoms = GetAllAtomsOfResidue(residue);
    
    for(PdbAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PdbAtom* atom = (*it);
        if(atom->GetAtomName().compare(atom_name) == 0)
            return atom;
    }
    return NULL;
}

PdbAtom* PdbFile::GetAtomOfResidueByAtomKey(string atom_key)
{
    vector<string> key_tokens = Split(atom_key, "_");
    int serial_number = gmml::ConvertString<int>(key_tokens.at(1));
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it = atom_cards.begin(); it != atom_cards.end(); it++)
    {
        PdbAtomCard* atom_card = (*it);
        PdbAtomCard::PdbAtomMap atom_map = atom_card->GetAtoms();
        if(atom_map[serial_number] != NULL)
            return atom_map[serial_number];
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it = heterogen_atom_cards.begin(); it != heterogen_atom_cards.end(); it++)
    {
        PdbHeterogenAtomCard* heterogen_atom_card = (*it);
        PdbHeterogenAtomCard::PdbHeterogenAtomMap heterogen_atom_map = heterogen_atom_card->GetHeterogenAtoms();
        if(heterogen_atom_map[serial_number] != NULL)
            return heterogen_atom_map[serial_number];
    }
    return NULL;
}

PdbAtom* PdbFile::GetAtomBySerialNumber(int serial_number)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModel* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
    for(PdbModelResidueSet::AtomCardVector::iterator it = atom_cards.begin(); it != atom_cards.end(); it++)
    {
        PdbAtomCard* atom_card = (*it);
        PdbAtomCard::PdbAtomMap atom_map = atom_card->GetAtoms();
        if(atom_map[serial_number] != NULL)
            return atom_map[serial_number];
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it = heterogen_atom_cards.begin(); it != heterogen_atom_cards.end(); it++)
    {
        PdbHeterogenAtomCard* heterogen_atom_card = (*it);
        PdbHeterogenAtomCard::PdbHeterogenAtomMap heterogen_atom_map = heterogen_atom_card->GetHeterogenAtoms();
        if(heterogen_atom_map[serial_number] != NULL)
            return heterogen_atom_map[serial_number];
    }
    return NULL;
}

vector<string> PdbFile::GetAllAtomNamesOfResidue(PdbResidue *residue, PdbFile::PdbResidueAtomsMap residue_atom_map)
{
    string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();
    PdbAtomVector* atoms = residue_atom_map[target_key];
    
    vector<string> atom_names;
    for(PdbAtomVector::iterator it = atoms->begin(); it != atoms->end(); it++)
    {
        PdbAtom* atom = (*it);
        atom_names.push_back(atom->GetAtomName());
    }
    return atom_names;
}

vector<string> PdbFile::GetAllAtomNamesOfResidue(PdbResidue *residue)
{
    PdbAtomVector atoms = GetAllAtomsOfResidue(residue);
    
    vector<string> atom_names;
    for(PdbAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PdbAtom* atom = (*it);
        atom_names.push_back(atom->GetAtomName());
    }
    return atom_names;
}
PdbFile::PdbPairVectorTerCardPositions PdbFile::GetAllTerCardPositions(vector<string> glycam_residue_names)
{
    vector<pair<char, int> > ter_card_positions = vector<pair<char, int> >();
    // After residues that has no tails or has more than or equal two tails
    PdbFile::PdbResidueVector residues = this->GetAllResiduesFromAtomCard();
    for(PdbFile::PdbResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        int dist = distance(residues.begin(), it);
        if(dist != residues.size() - 1)
        {
            PdbResidue* residue = (*it);
            string residue_name = residue->GetResidueName();
            if(find(glycam_residue_names.begin(), glycam_residue_names.end(), residue_name) != glycam_residue_names.end())
            {
                // No tail || has more than or equal two tails
                if(residue_name[0] == '0' || isalpha(residue_name[0]))
                {
                    ter_card_positions.push_back(make_pair(residue->GetResidueChainId(), residue->GetResidueSequenceNumber() + 1));
                }
            }
        }
    }
    return ter_card_positions;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbFile::SetPath(string pdb_path)
{
    path_ = pdb_path;
}
void PdbFile::SetHeader(PdbHeaderCard *header)
{
    header_ = new PdbHeaderCard();
    header_ = header;

}
void PdbFile::SetTitle(PdbTitleCard *title)
{
    title_ = new PdbTitleCard();
    title_ = title;
}
void PdbFile::SetCompound(PdbCompoundCard *compound)
{
    compound_ = new PdbCompoundCard();
    compound_ = compound;
}
void PdbFile::SetNumberOfModels(PdbNumModelCard *number_of_models)
{
    number_of_models_ = new PdbNumModelCard();
    number_of_models_ = number_of_models;
}
void PdbFile::SetModelType(PdbModelTypeCard *model_type)
{
    model_type_ = new PdbModelTypeCard();
    model_type_ = model_type;
}
void PdbFile::SetResiduesSequence(PdbResidueSequenceCard *residues_sequence)
{
    residues_sequence_ = new PdbResidueSequenceCard();
    residues_sequence_ = residues_sequence;
}
void PdbFile::SetResidueModification(PdbResidueModificationCard *residue_modification)
{
    residue_modification_ = new PdbResidueModificationCard();
    residue_modification_ = residue_modification;
}
void PdbFile::SetHeterogens(PdbHeterogenCard *heterogens)
{
    heterogens_ = new PdbHeterogenCard();
    heterogens_ = heterogens;
}
void PdbFile::SetHeterogensName(PdbHeterogenNameCard *heterogens_name)
{
    heterogens_name_ = new PdbHeterogenNameCard();
    heterogens_name_ = heterogens_name;
}
void PdbFile::SetHeterogenSynonyms(PdbHeterogenSynonymCard *heterogen_synonyms)
{
    heterogen_synonyms_ = new PdbHeterogenSynonymCard();
    heterogen_synonyms_ = heterogen_synonyms;
}
void PdbFile::SetFormulas(PdbFormulaCard *formulas)
{
    formulas_ = new PdbFormulaCard();
    formulas_ = formulas;
}
void PdbFile::SetHelixes(PdbHelixCard *helixes)
{
    helixes_ = new PdbHelixCard();
    helixes_ = helixes;
}
void PdbFile::SetSheets(PdbSheetCard *sheets)
{
    sheets_ = new PdbSheetCard();
    sheets_ = sheets;
}
void PdbFile::SetDisulfideBonds(PdbDisulfideBondCard *disulfide_bonds)
{
    disulfide_bonds_ = new PdbDisulfideBondCard();\
    disulfide_bonds_ = disulfide_bonds;
}
void PdbFile::SetLinks(PdbLinkCard *links)
{
    links_ = new PdbLinkCard();
    links_ = links;
}
void PdbFile::SetSites(PdbSiteCard *sites)
{
    sites_ = new PdbSiteCard();
    sites_ = sites;
}
void PdbFile::SetCrystallography(PdbCrystallographicCard *crystallography)
{
    crystallography_ = new PdbCrystallographicCard();
    crystallography_ = crystallography;
}
void PdbFile::SetOrigins(PdbOriginXnCard *origins)
{
    origins_ = new PdbOriginXnCard();
    origins_ = origins;
}
void PdbFile::SetScales(PdbScaleNCard *scales)
{
    scales_ = new PdbScaleNCard();
    scales_ = scales;
}
void PdbFile::SetMatrices(PdbMatrixNCard *matrices)
{
    matrices_ = new PdbMatrixNCard();
    matrices_ = matrices;
}
void PdbFile::SetModels(PdbModelCard *models)
{
    models_ = new PdbModelCard();
    models_ = models;
}
void PdbFile::SetConnectivities(PdbConnectCard *connectivities)
{
    connectivities_ = new PdbConnectCard();
    connectivities_ = connectivities;
}
void PdbFile::SetSerialNumberMapping(PdbSerialNumberMapping serial_number_mapping)
{
    serial_number_mapping_.clear();
    for(PdbSerialNumberMapping::iterator it = serial_number_mapping.begin(); it != serial_number_mapping.end(); it++)
    {
        int source_serial = (*it).first;
        int mapped_serial = (*it).second;
        serial_number_mapping_[source_serial] = mapped_serial;
    }
}
void PdbFile::SetSequenceNumberMapping(PdbSequenceNumberMapping sequence_number_mapping)
{
    sequence_number_mapping_.clear();
    for(PdbSequenceNumberMapping::iterator it = sequence_number_mapping.begin(); it != sequence_number_mapping.end(); it++)
    {
        int mapped_sequence = (*it).first;
        int source_sequence = (*it).second;
        sequence_number_mapping_[mapped_sequence] = source_sequence;
    }
}

void PdbFile::DeleteResidue(PdbResidue *residue)
{
    string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();
    
    PdbModelCard::PdbModelMap models = models_->GetModels();
    for(PdbModelCard::PdbModelMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModel* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms = PdbAtomCard::PdbAtomMap();
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();
            
            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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

            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
    map<string, PdbResidue* > pdb_residue_map = map<string, PdbResidue* >();
    for(PdbResidueVector::iterator it = target_residues.begin(); it != target_residues.end(); it++)
    {
        PdbResidue* residue = (*it);
        string target_residue_name = residue->GetResidueName();
        char target_residue_chain_id = residue->GetResidueChainId();
        int target_residue_sequence_number = residue->GetResidueSequenceNumber();
        char target_residue_insertion_code = residue->GetResidueInsertionCode();
        char target_residue_alternate_location = residue->GetResidueAlternateLocation();
        stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        string target_key = ss.str();
        pdb_residue_map[target_key] = residue;
    }

    PdbModelCard::PdbModelMap models = models_->GetModels();
    for(PdbModelCard::PdbModelMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModel* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms = PdbAtomCard::PdbAtomMap();
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
    string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();

    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModel* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms = PdbAtomCard::PdbAtomMap();
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
    map<string, PdbResidue* > pdb_residue_map = map<string, PdbResidue* >();
    for(PdbResidueVector::iterator it = target_residues.begin(); it != target_residues.end(); it++)
    {
        PdbResidue* residue = (*it);
        string target_residue_name = residue->GetResidueName();
        char target_residue_chain_id = residue->GetResidueChainId();
        int target_residue_sequence_number = residue->GetResidueSequenceNumber();
        char target_residue_insertion_code = residue->GetResidueInsertionCode();
        char target_residue_alternate_location = residue->GetResidueAlternateLocation();
        stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        string target_key = ss.str();
        pdb_residue_map[target_key] = residue;
    }

    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModel* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms = PdbAtomCard::PdbAtomMap();
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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

void PdbFile::DeleteAtom(PdbAtom* target_atom)
{
    string target_residue_name = target_atom->GetAtomResidueName();
    char target_residue_chain_id = target_atom->GetAtomChainId();
    int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
    char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
    char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();
    
    PdbModelCard::PdbModelMap models = models_->GetModels();
    for(PdbModelCard::PdbModelMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModel* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();
            
            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
                    string atom_name = atom->GetAtomName();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
                    string atom_name = atom->GetAtomName();
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

void PdbFile::DeleteAtoms(PdbAtomVector target_atoms)
{
    map<string, PdbAtom* > pdb_atom_map = map<string, PdbAtom* >();
    for(PdbAtomVector::iterator it = target_atoms.begin(); it != target_atoms.end(); it++)
    {
        PdbAtom* target_atom = (*it);
        string target_residue_name = target_atom->GetAtomResidueName();
        char target_residue_chain_id = target_atom->GetAtomChainId();
        int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
        char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
        char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
        stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        string target_key = ss.str();
        pdb_atom_map[target_key] = target_atom;
    }

    PdbModelCard::PdbModelMap models = models_->GetModels();
    for(PdbModelCard::PdbModelMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModel* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
                    PdbAtom* target_atom = pdb_atom_map[key];
                    string atom_name = atom->GetAtomName();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
                    PdbAtom* target_atom = pdb_atom_map[key];
                    string atom_name = atom->GetAtomName();
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

void PdbFile::DeleteAtomWithTheGivenModelNumber(PdbAtom* target_atom, int model_number)
{
    string target_residue_name = target_atom->GetAtomResidueName();
    char target_residue_chain_id = target_atom->GetAtomChainId();
    int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
    char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
    char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();

    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModel* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
                    string atom_name = atom->GetAtomName();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
                    string atom_name = atom->GetAtomName();
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

void PdbFile::DeleteAtomsWithTheGivenModelNumber(PdbAtomVector target_atoms, int model_number)
{
    map<string, PdbAtom* > pdb_atom_map = map<string, PdbAtom* >();
    for(PdbAtomVector::iterator it = target_atoms.begin(); it != target_atoms.end(); it++)
    {
        PdbAtom* target_atom = (*it);
        string target_residue_name = target_atom->GetAtomResidueName();
        char target_residue_chain_id = target_atom->GetAtomChainId();
        int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
        char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
        char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
        stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        string target_key = ss.str();
        pdb_atom_map[target_key] = target_atom;
    }

    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModel* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
                    PdbAtom* target_atom = pdb_atom_map[key];
                    string atom_name = atom->GetAtomName();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
                    PdbAtom* target_atom = pdb_atom_map[key];
                    string atom_name = atom->GetAtomName();
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

void PdbFile::UpdateResidueName(PdbResidue *residue, string updated_residue_name)
{
    string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();
    
    PdbModelCard::PdbModelMap models = models_->GetModels();
    for(PdbModelCard::PdbModelMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModel* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();
            
            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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

void PdbFile::UpdateResidueNameWithTheGivenModelNumber(PdbResidue *residue, string updated_residue_name, int model_number)
{
    string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    string target_key = ss.str();

    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModel* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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
            atom_card->SetAtoms(updated_atoms);
            atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
        }
        residue_set->SetAtoms(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = sss.str();
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

void PdbFile::InsertResidueBefore(PdbAtomCard* residue)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModelCard::PdbModelMap updated_models;
    for(PdbModelCard::PdbModelMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbModel* model = (*it).second;
        PdbModel* updated_model = new PdbModel();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard* updated_atom_card = new PdbAtomCard();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();
            bool located = true;

            PdbAtomCard::PdbAtomOrderVector ordered_atoms_of_residue = residue->GetOrderedAtoms();
            PdbAtom* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);

                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbAtom* previous_atom = (*(--it2));
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
                        // TODO: update coordinates with respect to it2
                        GeometryTopology::Coordinate::CoordinateVector coordinate_set = GeometryTopology::Coordinate::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                        base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, -1);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, -1);
                        base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, -1);
                        for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbAtom* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbAtom* new_atom = new PdbAtom(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        sequence_number_mapping_[sequence_number] = iNotSet;
                        sequence_number++;
                        offset++;
                        located = false;
                    }
                    if(!located)
                    {
                        PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
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

                        PdbAtom* previous_atom = (*(--it2));
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
                    PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
                    updated_atoms[serial_number] = updated_atom;
                    updated_atoms_vector.push_back(updated_atom);
                    serial_number++;
                    if(sequence_number != atom->GetAtomResidueSequenceNumber())
                        sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                }
            }
            updated_atom_card->SetAtoms(updated_atoms);
            updated_atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard* updated_heterogen_atom_card = new PdbHeterogenAtomCard();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* heterogen_atom = (*it2);
                PdbAtom* updated_heterogen_atom = new PdbAtom(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtoms(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
    }
    models_->SetModels(updated_models);
}

void PdbFile::InsertResidueBeforeWithTheGivenModelNumber(PdbAtomCard* residue, int model_number)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModel* model = models[model_number];
        PdbModelCard::PdbModelMap updated_models;
        PdbModel* updated_model = new PdbModel();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard* updated_atom_card = new PdbAtomCard();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();
            bool located = true;

            PdbAtomCard::PdbAtomOrderVector ordered_atoms_of_residue = residue->GetOrderedAtoms();
            PdbAtom* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);

                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbAtom* previous_atom = (*(--it2));
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
                        // TODO: update coordinates with respect to it2
                        GeometryTopology::Coordinate::CoordinateVector coordinate_set = GeometryTopology::Coordinate::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                        base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, -1);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, -1);
                        base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, -1);
                        for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbAtom* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbAtom* new_atom = new PdbAtom(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        sequence_number_mapping_[sequence_number] = iNotSet;
                        sequence_number++;
                        offset++;
                        located = false;
                    }
                    if(!located)
                    {
                        PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
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

                        PdbAtom* previous_atom = (*(--it2));
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
                    PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
                    updated_atoms[serial_number] = updated_atom;
                    updated_atoms_vector.push_back(updated_atom);
                    serial_number++;
                    if(sequence_number != atom->GetAtomResidueSequenceNumber())
                        sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                }
            }
            updated_atom_card->SetAtoms(updated_atoms);
            updated_atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard* updated_heterogen_atom_card = new PdbHeterogenAtomCard();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector  updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* heterogen_atom = (*it2);
                PdbAtom* updated_heterogen_atom = new PdbAtom(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtoms(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;

        models_->SetModels(updated_models);
    }
}

void PdbFile::InsertResidueAfter(PdbAtomCard* residue)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModelCard::PdbModelMap updated_models;
    for(PdbModelCard::PdbModelMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbModel* model = (*it).second;
        PdbModel* updated_model = new PdbModel();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard* updated_atom_card = new PdbAtomCard();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();
            bool located = false;
            
            PdbAtomCard::PdbAtomOrderVector ordered_atoms_of_residue = residue->GetOrderedAtoms();
            PdbAtom* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();
            
            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                
                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbAtom* previous_atom = (*(--it2));
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
                    PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
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
                        // TODO: update coordinates with respect to it2
                        GeometryTopology::Coordinate::CoordinateVector coordinate_set = GeometryTopology::Coordinate::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                        base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                        base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                        sequence_number++;
                        for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbAtom* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbAtom* new_atom = new PdbAtom(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        sequence_number_mapping_[sequence_number] = iNotSet;
                        sequence_number++;
                        offset++;
                        located = false;
                    }
                    if(!located)
                    {
                        if(it2 != ordered_atoms.begin())
                        {
                            PdbAtom* previous_atom = (*(--it2));
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
                        PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
                        updated_atoms[serial_number] = updated_atom;
                        updated_atoms_vector.push_back(updated_atom);
                        serial_number++;
                        if(sequence_number != atom->GetAtomResidueSequenceNumber())
                            sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    }
                }
                if(it2 == --ordered_atoms.end() && located)
                {
                    // TODO: update coordinates with respect to it2
                    GeometryTopology::Coordinate::CoordinateVector coordinate_set = GeometryTopology::Coordinate::CoordinateVector();
                    GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                    for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                    base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                    base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                    base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                    sequence_number++;
                    for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                    {
                        PdbAtom* atom_of_residue = (*it3);
                        int index = distance(ordered_atoms_of_residue.begin(), it3);
                        PdbAtom* new_atom = new PdbAtom(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                        atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                        atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                        atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                        atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge());
                        updated_atoms[serial_number] = new_atom;
                        updated_atoms_vector.push_back(new_atom);
                        serial_number++;
                    }
                    sequence_number_mapping_[sequence_number] = iNotSet;
                    located = false;
                    sequence_number++;
                    offset++;
                }
            }
            updated_atom_card->SetAtoms(updated_atoms);
            updated_atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard* updated_heterogen_atom_card = new PdbHeterogenAtomCard();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* heterogen_atom = (*it2);
                PdbAtom* updated_heterogen_atom = new PdbAtom(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtoms(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
    }
    models_->SetModels(updated_models);
}

void PdbFile::InsertResidueAfterWithTheGivenModelNumber(PdbAtomCard* residue, int model_number)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModel* model = models[model_number];
        PdbModelCard::PdbModelMap updated_models;
        PdbModel* updated_model = new PdbModel();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard* updated_atom_card = new PdbAtomCard();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap updated_atoms;
            PdbAtomCard::PdbAtomOrderVector updated_atoms_vector = PdbAtomCard::PdbAtomOrderVector();
            bool located = false;

            PdbAtomCard::PdbAtomOrderVector ordered_atoms_of_residue = residue->GetOrderedAtoms();
            PdbAtom* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);

                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbAtom* previous_atom = (*(--it2));
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
                    PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
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
                        // TODO: update coordinates with respect to it2
                        GeometryTopology::Coordinate::CoordinateVector coordinate_set = GeometryTopology::Coordinate::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                        base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                        base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                        sequence_number++;
                        for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbAtom* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbAtom* new_atom = new PdbAtom(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge());
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
                            PdbAtom* previous_atom = (*(--it2));
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
                        PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
                        updated_atoms[serial_number] = updated_atom;
                        updated_atoms_vector.push_back(updated_atom);
                        serial_number++;
                        if(sequence_number != atom->GetAtomResidueSequenceNumber())
                            sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    }
                }
                if(it2 == --ordered_atoms.end() && located)
                {
                    // TODO: update coordinates with respect to it2
                    GeometryTopology::Coordinate::CoordinateVector coordinate_set = GeometryTopology::Coordinate::CoordinateVector();
                    GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                    for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                    base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                    base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                    base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                    sequence_number++;
                    for(PdbAtomCard::PdbAtomOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                    {
                        PdbAtom* atom_of_residue = (*it3);
                        int index = distance(ordered_atoms_of_residue.begin(), it3);
                        PdbAtom* new_atom = new PdbAtom(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                        atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                        atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                        atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                        atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge());
                        updated_atoms[serial_number] = new_atom;
                        updated_atoms_vector.push_back(new_atom);
                        serial_number++;
                    }
                    sequence_number_mapping_[sequence_number] = iNotSet;
                    located = false;
                    sequence_number++;
                    offset++;
                }
            }
            updated_atom_card->SetAtoms(updated_atoms);
            updated_atom_card->SetOrderedAtoms(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard* updated_heterogen_atom_card = new PdbHeterogenAtomCard();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
            PdbHeterogenAtomCard::PdbHeterogenAtomMap updated_heterogen_atoms;
            PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbAtom* heterogen_atom = (*it2);
                PdbAtom* updated_heterogen_atom = new PdbAtom(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtoms(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
        models_->SetModels(updated_models);
    }
}
void PdbFile::SplitAtomCardOfModelCard(char split_point_chain_id, int split_point_sequence_number)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    PdbModelCard::PdbModelMap updated_models;
    for(PdbModelCard::PdbModelMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbModel* model = (*it).second;
        PdbModel* updated_model = new PdbModel();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard* updated_atom_card_first_part = new PdbAtomCard();
            PdbAtomCard* updated_atom_card_second_part = new PdbAtomCard();
            updated_atom_card_first_part->SetRecordName(atom_card->GetRecordName());
            updated_atom_card_second_part->SetRecordName(atom_card->GetRecordName());
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap atoms_first_part;
            PdbAtomCard::PdbAtomMap atoms_second_part;
            PdbAtomCard::PdbAtomOrderVector ordered_atoms_first_part = PdbAtomCard::PdbAtomOrderVector();
            PdbAtomCard::PdbAtomOrderVector ordered_atoms_second_part = PdbAtomCard::PdbAtomOrderVector();
            
            bool located = false;
            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                char residue_chain_id = atom->GetAtomChainId();
                int residue_sequence_number = atom->GetAtomResidueSequenceNumber();
                
                if(residue_chain_id == split_point_chain_id && residue_sequence_number == split_point_sequence_number)
                {
                    located = true;
                }
                if(located)
                {
                    serial_number++;
                    PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
                    atoms_second_part[serial_number] = updated_atom;
                    ordered_atoms_second_part.push_back(updated_atom);
                }
                if(!located)
                {
                    PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
                    atoms_first_part[serial_number] = updated_atom;
                    ordered_atoms_first_part.push_back(updated_atom);
                    serial_number++;
                }
            }
            updated_atom_card_first_part->SetAtoms(atoms_first_part);
            updated_atom_card_first_part->SetOrderedAtoms(ordered_atoms_first_part);
            updated_atom_card_second_part->SetAtoms(atoms_second_part);
            updated_atom_card_second_part->SetOrderedAtoms(ordered_atoms_second_part);
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
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        
        updated_residue_set->SetAtoms(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
    }
    models_->SetModels(updated_models);
}

void PdbFile::SplitAtomCardOfModelCardWithTheGivenModelNumber(char split_point_chain_id, int split_point_sequence_number, int model_number)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard::PdbModelMap updated_models;
        PdbModel* model = models[model_number];
        PdbModel* updated_model = new PdbModel();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            PdbAtomCard* updated_atom_card_first_part = new PdbAtomCard();
            PdbAtomCard* updated_atom_card_second_part = new PdbAtomCard();
            updated_atom_card_first_part->SetRecordName(atom_card->GetRecordName());
            updated_atom_card_second_part->SetRecordName(atom_card->GetRecordName());
            PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
            PdbAtomCard::PdbAtomMap atoms_first_part;
            PdbAtomCard::PdbAtomMap atoms_second_part;
            PdbAtomCard::PdbAtomOrderVector ordered_atoms_first_part = PdbAtomCard::PdbAtomOrderVector();
            PdbAtomCard::PdbAtomOrderVector ordered_atoms_second_part = PdbAtomCard::PdbAtomOrderVector();

            bool located = false;
            for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2);
                char residue_chain_id = atom->GetAtomChainId();
                int residue_sequence_number = atom->GetAtomResidueSequenceNumber();

                if(residue_chain_id == split_point_chain_id && residue_sequence_number == split_point_sequence_number)
                {
                    located = true;
                }
                if(located)
                {
                    serial_number++;
                    PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
                    atoms_second_part[serial_number] = updated_atom;
                    ordered_atoms_second_part.push_back(updated_atom);
                }
                if(!located)
                {
                    PdbAtom* updated_atom = new PdbAtom(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge());
                    atoms_first_part[serial_number] = updated_atom;
                    ordered_atoms_first_part.push_back(updated_atom);
                    serial_number++;
                }
            }
            updated_atom_card_first_part->SetAtoms(atoms_first_part);
            updated_atom_card_first_part->SetOrderedAtoms(ordered_atoms_first_part);
            updated_atom_card_second_part->SetAtoms(atoms_second_part);
            updated_atom_card_second_part->SetOrderedAtoms(ordered_atoms_second_part);
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
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();

        updated_residue_set->SetAtoms(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
        models_->SetModels(updated_models);
    }
}

void PdbFile::UpdateConnectCard()
{
    if(connectivities_ != NULL)
    {
        PdbConnectCard::BondedAtomsSerialNumbersMap bondedAtomsSerialNumbersMap = connectivities_->GetBondedAtomsSerialNumbers();
        PdbConnectCard::BondedAtomsSerialNumbersMap new_bonded_atoms_serial_numbers_map = PdbConnectCard::BondedAtomsSerialNumbersMap();
        for(PdbConnectCard::BondedAtomsSerialNumbersMap::iterator it = bondedAtomsSerialNumbersMap.begin(); it != bondedAtomsSerialNumbersMap.end(); it++)
        {
            int source_serial_number = (*it).first;
            vector<int> bonded_serial_numbers = (*it).second;
            if(serial_number_mapping_.find(source_serial_number) != serial_number_mapping_.end())
            {
                int new_source_serial_number = serial_number_mapping_[source_serial_number];
                new_bonded_atoms_serial_numbers_map[new_source_serial_number] = vector<int>();
                for(vector<int>::iterator it1 = bonded_serial_numbers.begin(); it1 != bonded_serial_numbers.end(); it1++)
                {
                    int bonded_serial_number = (*it1);
                    if(serial_number_mapping_.find(bonded_serial_number) != serial_number_mapping_.end())
                    {
                        int new_bonded_serial_number = serial_number_mapping_[bonded_serial_number];
                        new_bonded_atoms_serial_numbers_map[new_source_serial_number].push_back(new_bonded_serial_number);
                    }
                }
            }
        }
        connectivities_->SetBondedAtomsSerialNumbers(new_bonded_atoms_serial_numbers_map);
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
bool PdbFile::Read(ifstream &in_file)
{
    if(!this->ParseCards(in_file))
        return false;
}

bool PdbFile::ParseCards(ifstream &in_stream)
{
    string line;
    
    /// Unable to read file
    if (!getline(in_stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
        cout << "Wrong input file format" << endl;
        throw PdbFileProcessingException("Error reading file");
    }

    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("HEADER") == 0)
    {
        if(!ParseHeaderCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("OBSLTE") == 0)
    {
        if(!ParseObsoleteCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("TITLE") == 0)
    {
        if(!ParseTitleCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("SPLIT") == 0)
    {
        if(!ParseSplitCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("CAVEAT") == 0)
    {
        if(!ParseCaveatCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("COMPND") == 0)
    {
        if(!ParseCompoundCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("SOURCE") == 0)
    {
        if(!ParseSourceCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("KEYWDS") == 0)
    {
        if(!ParseKeywordCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("EXPDTA") == 0)
    {
        if(!ParseExpirationDateCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("NUMMDL") == 0)
    {
        if(!ParseNumModelCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("MDLTYP") == 0)
    {
        if(!ParseModelTypeCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("AUTHOR") == 0)
    {
        if(!ParseAuthorCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("REVDAT") == 0)
    {
        if(!ParseRevisionDateCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("SPRSDE") == 0)
    {
        if(!ParseSupersededEntriesCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("JRNL") == 0)
    {
        if(!ParseJournalCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("REMARK") == 0)
    {
        if(!ParseRemarkCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("DBREF") != string::npos)
    {
        if(!ParseDatabaseReferenceCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("SEQADV") == 0)
    {
        if(!ParseSequenceAdvancedCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("SEQRES") == 0)
    {
        if(!ParseSequenceResidueCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("MODRES") == 0)
    {
        if(!ParseModificationResidueCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("HET") == 0)
    {
        if(!ParseHeterogenCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("HETNAM") == 0)
    {
        if(!ParseHeterogenNameCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("HETSYN") == 0)
    {
        if(!ParseHeterogenSynonymCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("FORMUL") == 0)
    {
        if(!ParseFormulaCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("HELIX") == 0)
    {
        if(!ParseHelixCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("SHEET") == 0)
    {
        if(!ParseSheetCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("SSBOND") == 0)
    {
        if(!ParseDisulfideBondCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("LINK") == 0)
    {
        if(!ParseLinkCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("CISPEP") == 0)
    {
        if(!ParseCISPeptideCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("SITE") == 0)
    {
        if(!ParseSiteCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("CRYST1") == 0)
    {
        if(!ParseCrystallographyCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("ORIGX") != string::npos)
    {
        if(!ParseOriginCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("SCALE") != string::npos)
    {
        if(!ParseScaleCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("MTRIX") != string::npos)
    {
        if(!ParseMatrixCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("MODEL") == 0 || record_name.compare("ATOM") == 0 || record_name.compare("HETATM") == 0)
    {
        if(!ParseModelCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.compare("CONECT") == 0)
    {
        if(!ParseConnectivityCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);

    if(record_name.compare("MODEL") == 0)
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Multiple connect card between model cards --> Unexpected entry");
        cout << "Multiple connect card between model cards --> Unexpected entry" << endl;
        return false;
    }
    if(record_name.compare("MASTER") == 0)
    {
        if(!ParseMasterCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("END") != string::npos || record_name.compare("END") == 0)
    {
        if(!ParseEndCard(in_stream, line))
            return false;
        return true;
    }
    else
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
        cout << "Wrong input file format" << endl;
        stringstream ss;
        ss << record_name << " is an Unknown record name.";
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        cout << ss.str() << endl;
        return false;
    }
    return true;
}

bool PdbFile::ParseHeaderCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Header card corupption");
        cout << "Header card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("HEADER") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Header card corruption");
            cout << "Header card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    header_ = new PdbHeaderCard(stream_block);
    return true;
}

bool PdbFile::ParseObsoleteCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Obsolete card corruption");
        cout << "Obsolete card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("OBSLTE") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Obsolete card corruption");
            cout << "Obsolete card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseTitleCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Title card corruption");
        cout << "Title card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("TITLE") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Title card corruption");
            cout << "Title card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    title_ = new PdbTitleCard(stream_block);
    return true;
}

bool PdbFile::ParseSplitCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Split card corruption" );
        cout << "Split card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("SPLIT") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Split card corruption" );
            cout << "Split card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseCaveatCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Caveat card corruption" );
        cout << "Caveat card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("CAVEAT") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Caveat card corruption" );
            cout << "Caveat card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseCompoundCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Compound card corruption" );
        cout << "Compound card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("COMPND") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Compound card corruption" );
            cout << "Compound card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    compound_ = new PdbCompoundCard(stream_block);
    return true;
}

bool PdbFile::ParseSourceCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Source card corruption" );
        cout << "Source card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("SOURCE") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Source card corruption" );
            cout << "Source card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseKeywordCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Keyword card corruption" );
        cout << "Keyword card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("KEYWDS") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Keyword card corruption" );
            cout << "Keyword card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseExpirationDateCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Expiration date card corruption" );
        cout << "Expiration date card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("EXPDTA") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Expiration date card corruption" );
            cout << "Expiration date card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseNumModelCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Number of model card corruption" );
        cout << "Number of model card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("NUMMDL") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Number of model card corruption" );
            cout << "Number of model card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    number_of_models_ = new PdbNumModelCard(stream_block);
    return true;
}

bool PdbFile::ParseModelTypeCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Model type card corruption" );
        cout << "Model type card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("MDLTYP") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Model type card corruption" );
            cout << "Model type card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    model_type_ = new PdbModelTypeCard(stream_block);
    return true;
}

bool PdbFile::ParseAuthorCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Author card corruption" );
        cout << "Author card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("AUTHOR") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Author card corruption" );
            cout << "Author card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseRevisionDateCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Revision date card corruption" );
        cout << "Revision date card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("REVDAT") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Revision date card corruption" );
            cout << "Revision date card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseSupersededEntriesCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Superseded entries card corruption" );
        cout << "Superseded entries card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("SPRSDE") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Superseded entries card corruption" );
            cout << "Superseded entries card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseJournalCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Journal card corruption" );
        cout << "Journal card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("JRNL") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Journal card corruption" );
            cout << "Journal card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseRemarkCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Remark card corruption" );
        cout << "Remark card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("REMARK") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Remark card corruption" );
            cout << "Remark card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseDatabaseReferenceCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "database reference card corruption" );
        cout << "database reference card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("DBREF") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "database reference card corruption" );
            cout << "database reference card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseSequenceAdvancedCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sequence advanced card corruption" );
        cout << "Sequence advanced card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("SEQADV") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sequence advanced card corruption" );
            cout << "Sequence advanced card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseSequenceResidueCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sequence residue card corruption" );
        cout << "Sequence residue card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("SEQRES") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sequence residue card corruption" );
            cout << "Sequence residue card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    residues_sequence_ = new PdbResidueSequenceCard(stream_block);
    return true;
}

bool PdbFile::ParseModificationResidueCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Modification residue card corruption" );
        cout << "Modification residue card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("MODRES") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Modification residue card corruption" );
            cout << "Modification residue card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    residue_modification_ = new PdbResidueModificationCard(stream_block);
    return true;
}

bool PdbFile::ParseHeterogenCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen card corruption" );
        cout << "Heterogen card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("HET") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen card corruption" );
            cout << "Heterogen card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    heterogens_ = new PdbHeterogenCard(stream_block);
    return true;
}

bool PdbFile::ParseHeterogenNameCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen name card corruption" );
        cout << "Heterogen name card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("HETNAM") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen name card corruption" );
            cout << "Heterogen name card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    heterogens_name_ = new PdbHeterogenNameCard(stream_block);
    return true;
}

bool PdbFile::ParseHeterogenSynonymCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen synonym card corruption" );
        cout << "Heterogen synonym card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("HETSYN") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen synonym card corruption" );
            cout << "Heterogen synonym card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    heterogen_synonyms_ = new PdbHeterogenSynonymCard(stream_block);
    return true;
}

bool PdbFile::ParseFormulaCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Formula card corruption" );
        cout << "Formula card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("FORMUL") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Formula card corruption" );
            cout << "Formula card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    formulas_ = new PdbFormulaCard(stream_block);
    return true;
}

bool PdbFile::ParseHelixCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Helix card corruption" );
        cout << "Helix card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("HELIX") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Helix card corruption" );
            cout << "Helix card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    helixes_ = new PdbHelixCard(stream_block);
    return true;
}

bool PdbFile::ParseSheetCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sheet card corruption" );
        cout << "Sheet card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name.compare("SHEET") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sheet card corruption" );
            cout << "Sheet card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    sheets_ = new PdbSheetCard(stream_block);
    return true;
}

bool PdbFile::ParseDisulfideBondCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Disulfide bond card corruption" );
        cout << "Disulfide bond card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("SSBOND") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Disulfide bond card corruption" );
            cout << "Disulfide bond card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    disulfide_bonds_ = new PdbDisulfideBondCard(stream_block);
    return true;
}

bool PdbFile::ParseLinkCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Link card corruption" );
        cout << "Link card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("LINK") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Link card corruption" );
            cout << "Link card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    links_ = new PdbLinkCard(stream_block);
    return true;
}

bool PdbFile::ParseCISPeptideCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "CIS peptide card corruption" );
        cout << "CIS peptide card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("CISPEP") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "CIS peptide card corruption" );
            cout << "CIS peptide card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseSiteCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Site card corruption" );
        cout << "Site card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("SITE") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Site card corruption" );
            cout << "Site card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    sites_ = new PdbSiteCard(stream_block);
    return true;
}

bool PdbFile::ParseCrystallographyCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Crystallography card corruption" );
        cout << "Crystallography card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("CRYST1") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Crystallography card corruption" );
            cout << "Crystallography card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    crystallography_ = new PdbCrystallographicCard(stream_block);
    return true;
}

bool PdbFile::ParseOriginCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Origin card corruption" );
        cout << "Origin card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,5);
    record_name = Trim(record_name);
    
    while(record_name.compare("ORIGX") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,5);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Origin card corruption" );
            cout << "Origin card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    origins_ = new PdbOriginXnCard(stream_block);
    return true;
}

bool PdbFile::ParseScaleCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Scale card corruption" );
        cout << "Scale card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,5);
    record_name = Trim(record_name);
    
    while(record_name.compare("SCALE") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,5);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Scale card corruption" );
            cout << "Scale card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    scales_ = new PdbScaleNCard(stream_block);
    return true;
}

bool PdbFile::ParseMatrixCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Matrix card corruption" );
        cout << "Matrix card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,5);
    record_name = Trim(record_name);
    
    while(record_name.compare("MTRIX") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,5);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Matrix card corruption" );
            cout << "Matrix card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    matrices_ = new PdbMatrixNCard(stream_block);
    return true;
}

bool PdbFile::ParseModelCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Model card corruption" );
        cout << "Model card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name.compare("MODEL") == 0 || record_name.compare("ATOM") == 0 || record_name.compare("ANISOU") == 0
          || record_name.compare("TER") == 0 || record_name.compare("HETATM") == 0 || record_name.compare("ENDMDL") == 0)
        //          || record_name.find("TER") != string::npos || record_name.find("ENDMDL") != string::npos)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Model card corruption" );
            cout << "Model card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    // Model card
    //    gmml::log(__LINE__, __FILE__,  gmml::ERR, stream_block.str();
    models_ = new PdbModelCard(stream_block);
    return true;
}

bool PdbFile::ParseConnectivityCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Connectivity card corruption" );
        cout << "Connectivity card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("CONECT") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Connectivity card corruption" );
            cout << "Connectivity card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    
    connectivities_ = new PdbConnectCard(stream_block);
    return true;
}

bool PdbFile::ParseMasterCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Master card corruption" );
        cout << "Master card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
        cout << "Wrong input file format" << endl;
        return false;
    }
    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    
    while(record_name.compare("MASTER") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Master card corruption" );
            cout << "Master card corruption" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
            cout << "Wrong input file format" << endl;
            return false;
        }
    }
    return true;
}

bool PdbFile::ParseEndCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "End of file" );
        cout << "End of file" << endl;
        return true;
    }

    line = ExpandLine(line, iPdbLineLength);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    while(record_name.find("END") != string::npos || record_name.compare("END") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            line = ExpandLine(line, iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::INF, "End of file" );
            return true;
        }
    }
    return true;
}

void PdbFile::Write(const std::string& pdb_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(pdb_file.c_str());
    }
    catch(...)
    {
        throw PdbFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->ResolveCards(out_file);
    }
    catch(...)
    {
        out_file.close();            /// Close the parameter files
    }
}

void PdbFile::WriteWithTheGivenModelNumber(const std::string& pdb_file, int model_number)
{
    std::ofstream out_file;
    try
    {
        out_file.open(pdb_file.c_str());
    }
    catch(...)
    {
        throw PdbFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->ResolveCardsWithTheGivenModelNumber(out_file, model_number);
    }
    catch(...)
    {
        out_file.close();            /// Close the parameter files
    }
}

void PdbFile::ResolveCards(std::ofstream& out_stream)
{
    if(this->header_ != NULL)
    {
        this->ResolveHeaderCard(out_stream);
    }
    if(this->title_ != NULL)
    {
        this->ResolveTitleCard(out_stream);
    }
    if(this->compound_ != NULL)
    {
        this->ResolveCompoundCard(out_stream);
    }
    if(this->number_of_models_ != NULL)
    {
        this->ResolveCompoundCard(out_stream);
    }
    if(this->model_type_ != NULL)
    {
        this->ResolveNumModelCard(out_stream);
    }
    if(this->residues_sequence_ != NULL)
    {
        this->ResolveSequenceResidueCard(out_stream);
    }
    if(this->residue_modification_ != NULL)
    {
        this->ResolveModificationResidueCard(out_stream);
    }
    if(this->heterogens_ != NULL)
    {
        this->ResolveHeterogenCard(out_stream);
    }
    if(this->heterogens_name_ != NULL)
    {
        this->ResolveHeterogenNameCard(out_stream);
    }
    if(this->heterogen_synonyms_ != NULL)
    {
        this->ResolveHeterogenSynonymCard(out_stream);
    }
    if(this->formulas_ != NULL)
    {
        this->ResolveFormulaCard(out_stream);
    }
    if(this->helixes_ != NULL)
    {
        this->ResolveHelixCard(out_stream);
    }
    if(this->sheets_ != NULL)
    {
        this->ResolveSheetCard(out_stream);
    }
    if(this->disulfide_bonds_ != NULL)
    {
        this->ResolveDisulfideBondCard(out_stream);
    }
    if(this->links_ != NULL)
    {
        this->ResolveLinkCard(out_stream);
    }
    if(this->sites_ != NULL)
    {
        this->ResolveSiteCard(out_stream);
    }
    if(this->crystallography_ != NULL)
    {
        this->ResolveCrystallographyCard(out_stream);
    }
    if(this->origins_ != NULL)
    {
        this->ResolveOriginCard(out_stream);
    }
    if(this->scales_ != NULL)
    {
        this->ResolveScaleCard(out_stream);
    }
    if(this->matrices_ != NULL)
    {
        this->ResolveMatrixCard(out_stream);
    }
    if(this->models_ != NULL)
    {
        this->ResolveModelCard(out_stream);
    }
    if(this->connectivities_ != NULL)
    {
        this->ResolveConnectivityCard(out_stream);
    }
    this->ResolveEndCard(out_stream);
}

void PdbFile::ResolveCardsWithTheGivenModelNumber(std::ofstream& out_stream, int model_number)
{
    if(this->header_ != NULL)
    {
        this->ResolveHeaderCard(out_stream);
    }
    if(this->title_ != NULL)
    {
        this->ResolveTitleCard(out_stream);
    }
    if(this->compound_ != NULL)
    {
        this->ResolveCompoundCard(out_stream);
    }
    if(this->number_of_models_ != NULL)
    {
        this->ResolveCompoundCard(out_stream);
    }
    if(this->model_type_ != NULL)
    {
        this->ResolveNumModelCard(out_stream);
    }
    if(this->residues_sequence_ != NULL)
    {
        this->ResolveSequenceResidueCard(out_stream);
    }
    if(this->residue_modification_ != NULL)
    {
        this->ResolveModificationResidueCard(out_stream);
    }
    if(this->heterogens_ != NULL)
    {
        this->ResolveHeterogenCard(out_stream);
    }
    if(this->heterogens_name_ != NULL)
    {
        this->ResolveHeterogenNameCard(out_stream);
    }
    if(this->heterogen_synonyms_ != NULL)
    {
        this->ResolveHeterogenSynonymCard(out_stream);
    }
    if(this->formulas_ != NULL)
    {
        this->ResolveFormulaCard(out_stream);
    }
    if(this->helixes_ != NULL)
    {
        this->ResolveHelixCard(out_stream);
    }
    if(this->sheets_ != NULL)
    {
        this->ResolveSheetCard(out_stream);
    }
    if(this->disulfide_bonds_ != NULL)
    {
        this->ResolveDisulfideBondCard(out_stream);
    }
    if(this->links_ != NULL)
    {
        this->ResolveLinkCard(out_stream);
    }
    if(this->sites_ != NULL)
    {
        this->ResolveSiteCard(out_stream);
    }
    if(this->crystallography_ != NULL)
    {
        this->ResolveCrystallographyCard(out_stream);
    }
    if(this->origins_ != NULL)
    {
        this->ResolveOriginCard(out_stream);
    }
    if(this->scales_ != NULL)
    {
        this->ResolveScaleCard(out_stream);
    }
    if(this->matrices_ != NULL)
    {
        this->ResolveMatrixCard(out_stream);
    }
    if(this->models_ != NULL)
    {
        this->ResolveModelCardWithTheGivenModelNumber(out_stream, model_number);
    }
    if(this->connectivities_ != NULL)
    {
        this->ResolveConnectivityCard(out_stream);
    }
    this->ResolveEndCard(out_stream);
}

void PdbFile::ResolveHeaderCard(std::ofstream& stream)
{
    stream << left << setw(6) << header_->GetRecordName()
           << left << setw(4) << " "
           << left << setw(40) << header_->GetClassification()
           << left << setw(9) << header_->GetDepositionDate()
           << left << setw(3) << " "
           << right << setw(4) << header_->GetIdentifierCode()
           << left << setw(14) << " "
           << endl;
}

void PdbFile::ResolveObsoleteCard(std::ofstream& stream)
{
}

void PdbFile::ResolveTitleCard(std::ofstream& stream)
{
    const int MAX_TITLE_LENGTH_IN_LINE = 70;
    stream << left << setw(6) << title_->GetRecordName()
           << left << setw(2) << " ";
    if((int)title_->GetTitle().length() > MAX_TITLE_LENGTH_IN_LINE)
    {
        stream << right << setw(2) << " "
               << left << setw(70) << title_->GetTitle().substr(0,MAX_TITLE_LENGTH_IN_LINE)
               << endl;
        
        int counter = ceil((double)(title_->GetTitle().length()) / MAX_TITLE_LENGTH_IN_LINE);
        for(int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << left << setw(6) << title_->GetRecordName()
                       << left << setw(2) << " "
                       << right << setw(2) << i
                       << left << setw(70) << title_->GetTitle().substr(MAX_TITLE_LENGTH_IN_LINE*(i-1), MAX_TITLE_LENGTH_IN_LINE)
                       << endl;
            }
            else
            {
                stream << left << setw(6) << title_->GetRecordName()
                       << left << setw(2) << " "
                       << right << setw(2) << i
                       << left << setw(70) << title_->GetTitle().substr(MAX_TITLE_LENGTH_IN_LINE*(i-1), title_->GetTitle().length()-MAX_TITLE_LENGTH_IN_LINE*(i-1))
                       << endl;
            }
        }
    }
    else
    {
        stream << right << setw(2) << " "
               << left << setw(70) << title_->GetTitle()
               << endl;
    }
}

void PdbFile::ResolveSplitCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveCaveatCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveCompoundCard(std::ofstream& stream)
{
    const int MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE = 70;
    stream << left << setw(6) << compound_->GetRecordName()
           << left << setw(1) << " "
           << right << setw(3) << " ";
    
    PdbCompoundCard::PdbCompoundSpecificationMap compound_specification_map = compound_->GetCompoundSpecifications();
    if((*(compound_specification_map.begin())).second->GetMoleculeId() != "")
    {
        stringstream ss;
        ss << "MOL_ID: " << (*(compound_specification_map.begin())).second->GetMoleculeId() << ";";
        stream << left << setw(70) << ss.str() << endl;
    }
    else
    {
        stringstream ss;
        ss << " UNKNOWN;";
        stream << left << setw(70) << ss.str() << endl;
    }
    bool first = true;
    int counter = 2;
    for(PdbCompoundCard::PdbCompoundSpecificationMap::iterator it = compound_specification_map.begin(); it != compound_specification_map.end(); it++)
    {
        PdbCompoundSpecification* compound_specification = (*it).second;
        if(!first)
        {
            if(compound_specification->GetMoleculeId() != "")
            {
                stringstream ss;
                ss << " MOL_ID: " << compound_specification->GetMoleculeId() << ";";
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
            else
            {
                stringstream ss;
                ss << " UNKNOWN;";
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
        }
        
        /// Molecule name specification
        if(compound_specification->GetMoleculeName() != "")
        {
            stringstream molecule_name;
            molecule_name << " MOLECULE: " << compound_specification->GetMoleculeName() << ";";
            int length = molecule_name.str().length();
            
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                stringstream ss;
                ss << molecule_name.str();
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        stringstream ss;
                        ss << molecule_name.str().substr((i-1)*(MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE), MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str()
                               << endl;
                        counter++;
                    }
                    else
                    {
                        stringstream ss;
                        ss << molecule_name.str().substr((i-1)*(MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE), length - (i-1)*(MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE));
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str()
                               << endl;
                        counter++;
                    }
                }
                
            }
        }
        
        /// Molecule chain ids specification
        if(compound_specification->GetChainIds().size() > 0)
        {
            vector<string> chain_ids = compound_specification->GetChainIds();
            stringstream chain_id;
            chain_id << " CHAIN: ";
            for(vector<string>::iterator it1 = chain_ids.begin(); it1 != chain_ids.end(); it1++)
            {
                if(it1 < chain_ids.end()-1)
                    chain_id << (*it1) << ",";
                else
                    chain_id << (*it1);
            }
            chain_id << ";";
            int length = chain_id.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                stringstream ss;
                ss << chain_id.str();
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        stringstream ss;
                        ss << chain_id.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                    else
                    {
                        stringstream ss;
                        ss << chain_id.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                }
            }
        }
        
        /// Fragment specification
        if(compound_specification->GetFragment() != "")
        {
            stringstream fragment;
            fragment << " FRAGMENT: " << compound_specification->GetFragment() << ";";
            int length = fragment.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                stringstream ss;
                ss << ss.str();
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        stringstream ss;
                        ss << ss.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                    else
                    {
                        stringstream ss;
                        ss << ss.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                }
            }
        }
        
        /// Molecule synonyms specification
        if(compound_specification->GetMoleculeSynonyms().size() > 0)
        {
            vector<string> molecule_synonyms = compound_specification->GetMoleculeSynonyms();
            stringstream synonyms;
            synonyms << " SYNONYM: ";
            for(vector<string>::iterator it1 = molecule_synonyms.begin(); it1 != molecule_synonyms.end(); it1++)
            {
                if(it1 < molecule_synonyms.end()-1)
                    synonyms << (*it1) << ",";
                else
                    synonyms << (*it1);
            }
            synonyms << ";";
            int length = synonyms.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                stringstream ss;
                ss << synonyms.str();
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        stringstream ss;
                        ss << synonyms.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                    else
                    {
                        stringstream ss;
                        ss << synonyms.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                }
            }
        }
        
        /// Enzyme commission numbers specification
        if(compound_specification->GetEnzymeCommissionNumbers().size() > 0)
        {
            vector<string> enzyme_commission_numbers = compound_specification->GetEnzymeCommissionNumbers();
            stringstream commission_numbers;
            commission_numbers << " EC: ";
            for(vector<string>::iterator it1 = enzyme_commission_numbers.begin(); it1 != enzyme_commission_numbers.end(); it1++)
            {
                if(it1 < enzyme_commission_numbers.end()-1)
                    commission_numbers << (*it1) << ",";
                else
                    commission_numbers << (*it1);
            }
            commission_numbers << ";";
            int length = commission_numbers.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                stringstream ss;
                ss << commission_numbers.str();
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        stringstream ss;
                        ss << commission_numbers.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                    else
                    {
                        stringstream ss;
                        ss << commission_numbers.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                }
            }
        }
        
        /// Engineered specification
        if(compound_specification->GetIsEngineered())
        {
            stringstream ss;
            ss << " ENGINEERED: YES;";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }
        
        /// Mutation specification
        if(compound_specification->GetHasMutation())
        {
            stringstream ss;
            ss << " MUTATION: YES;";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }
        
        /// Other comments specification
        if(compound_specification->GetComments() != "")
        {
            stringstream comments;
            comments << " OTHER_DETAILS: " << compound_specification->GetComments() << ";";
            int length = comments.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                stringstream ss;
                ss << comments.str();
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        stringstream ss;
                        ss << comments.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                    else
                    {
                        stringstream ss;
                        ss << comments.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << left << setw(6) << compound_->GetRecordName()
                               << left << setw(1) << " "
                               << right << setw(3) << counter
                               << left << setw(70) << ss.str() << endl;
                        counter++;
                    }
                }
            }
        }
        
        first = false;
        
    }
    
}

void PdbFile::ResolveSourceCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveKeywordCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveExpirationDateCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveNumModelCard(std::ofstream& stream)
{
    stream << left << setw(6) << number_of_models_->GetRecordName()
           << left << setw(4) << " ";
    if(number_of_models_->GetNumberOfModels() != iNotSet)
        stream << right << setw(4) << number_of_models_->GetNumberOfModels();
    else
        stream << right << setw(4) << " ";
    stream << left << setw(66) << " "
           << endl;
}

void PdbFile::ResolveModelTypeCard(std::ofstream& stream)
{
    stream << left << setw(6) << model_type_->GetRecordName()
           << left << setw(2) << " ";
    stringstream ss;
    for(vector<string>::iterator it = model_type_->GetComments().begin(); it != model_type_->GetComments().end(); it++)
    {
        if(it != model_type_->GetComments().end() - 1)
        {
            ss << (*it) << "; ";
        }
        else
        {
            ss << (*it);
        }
    }
    if(ss.str().length() > 70)
    {
        stream << right << setw(2) << " "
               << left << setw(70) << ss.str().substr(0,70)
               << endl;
        
        int counter = ceil((double)(ss.str().length()) / 70);
        for(int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << left << setw(6) << model_type_->GetRecordName()
                       << left << setw(2) << " "
                       << right << setw(2) << i
                       << left << setw(70) << ss.str().substr(70*(i-1), 70)
                       << endl;
            }
            else
            {
                stream << left << setw(6) << model_type_->GetRecordName()
                       << left << setw(2) << " "
                       << right << setw(2) << i
                       << left << setw(70) << ss.str().substr(70*(i-1), ss.str().length() - (i-1)*70)
                       << endl;
            }
        }
    }
    else
    {
        stream << right << setw(2) << " "
               << left << setw(70) << ss.str()
               << endl;
    }
}

void PdbFile::ResolveAuthorCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveRevisionDateCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveSupersededEntriesCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveJournalCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveRemarkCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveDatabaseReferenceCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveSequenceAdvancedCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveSequenceResidueCard(std::ofstream& stream)
{
    PdbResidueSequenceCard::ResidueSequenceMap residue_sequence_map = residues_sequence_->GetResidueSequenceChain();
    for(PdbResidueSequenceCard::ResidueSequenceMap::iterator it = residue_sequence_map.begin(); it != residue_sequence_map.end(); it++)
    {
        
        PdbResidueSequence* residue_sequence = (*it).second;
        int serial_number = 1;
        const int MAX_RESIDUE_IN_SINGLE_LINE = 13;
        vector<string> residue_names = residue_sequence->GetResidueNames();
        if(residue_sequence->GetNumberOfResidues() <= MAX_RESIDUE_IN_SINGLE_LINE)
        {
            stringstream ss;
            for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
            {
                ss << right << setw(1) << " "
                   << right << setw(3) << (*it1);
            }
            for(int i = 0; i < MAX_RESIDUE_IN_SINGLE_LINE - residue_sequence->GetNumberOfResidues(); i++)
            {
                ss << right << setw(1) << " "
                   << right << setw(3) << " ";
            }
            stream << left << setw(6) << residues_sequence_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << serial_number
                   << right << setw(1) << " "
                   << right << setw(1) << residue_sequence->GetChainId()
                   << right << setw(1) << " ";
            if(residue_sequence->GetNumberOfResidues() != iNotSet)
                stream << right << setw(4) << residue_sequence->GetNumberOfResidues();
            else
                stream << right << setw(4) << " ";
            stream << right << setw(1) << " "
                   << right << setw(52) << ss.str()
                   << right << setw(10) << " "
                   << endl;
        }
        else
        {
            int number_of_lines = ceil((double)(residue_sequence->GetNumberOfResidues()) / MAX_RESIDUE_IN_SINGLE_LINE);
            for(int i = 0; i < number_of_lines; i++)
            {
                stringstream ss;
                if(i != number_of_lines - 1)
                {
                    for(vector<string>::iterator it1 = residue_names.begin() + i * MAX_RESIDUE_IN_SINGLE_LINE;
                        it1 != residue_names.begin() + (i+1) * MAX_RESIDUE_IN_SINGLE_LINE; it1++)
                    {
                        ss << right << setw(1) << " "
                           << right << setw(3) << (*it1);
                    }
                }
                else
                {
                    for(vector<string>::iterator it1 = residue_names.begin() + i * MAX_RESIDUE_IN_SINGLE_LINE;
                        it1 != residue_names.end(); it1++)
                    {
                        ss << right << setw(1) << " "
                           << right << setw(3) << (*it1);
                    }
                    for(int i = 0; i < MAX_RESIDUE_IN_SINGLE_LINE - residue_sequence->GetNumberOfResidues() % MAX_RESIDUE_IN_SINGLE_LINE; i++)
                    {
                        ss << right << setw(1) << " "
                           << right << setw(3) << " ";
                    }
                }
                stream << left << setw(6) << residues_sequence_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << serial_number
                       << right << setw(1) << " "
                       << right << setw(1) << residue_sequence->GetChainId()
                       << right << setw(1) << " ";
                if(residue_sequence->GetNumberOfResidues() != iNotSet)
                    stream << right << setw(4) << residue_sequence->GetNumberOfResidues();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << " "
                       << right << setw(52) << ss.str()
                       << right << setw(10) << " "
                       << endl;
                serial_number++;
            }
        }
        
    }
}

void PdbFile::ResolveModificationResidueCard(std::ofstream& stream)
{
    PdbResidueModificationCard::ResidueModificationMap residue_modification_map = residue_modification_->GetResidueModifications();
    for(PdbResidueModificationCard::ResidueModificationMap::iterator it = residue_modification_map.begin(); it != residue_modification_map.end(); it++)
    {
        PdbResidueModification* residue_modification = (*it).second;
        stream << left << setw(6) << residue_modification_->GetRecordName()
               << left << setw(1) << " "
               << right << setw(4) << residue_modification->GetIdCode()
               << left << setw(1) << " "
               << right << setw(3) << residue_modification->GetResidueName()
               << left << setw(1) << " "
               << right << setw(1) << residue_modification->GetChainId()
               << left << setw(1) << " ";
        if(residue_modification->GetSequenceNumber() != iNotSet)
            stream << right << setw(4) << residue_modification->GetSequenceNumber();
        else
            stream << right << setw(4) << " ";
        stream << right << setw(1) << residue_modification->GetInsertionCode()
               << left << setw(1) << " "
               << right << setw(3) << residue_modification->GetStandardResidueName()
               << left << setw(2) << " "
               << left << setw(41) << residue_modification->GetDscr()
               << left << setw(10) << " "
               << endl;
    }
}

void PdbFile::ResolveHeterogenCard(std::ofstream& stream)
{
    PdbHeterogenCard::HeterogenMap heterogen_map = heterogens_->GetHeterogens();
    for(PdbHeterogenCard::HeterogenMap::iterator it = heterogen_map.begin(); it != heterogen_map.end(); it++)
    {
        
        PdbHeterogen* heterogen = (*it).second;
        stream << left << setw(6) << heterogens_->GetRecordName()
               << left << setw(1) << " "
               << right << setw(3) << heterogen->GetHeterogenId()
               << left << setw(2) << " "
               << right << setw(1) << heterogen->GetChainId();
        if(heterogen->GetSequenceNumber() != iNotSet)
            stream << right << setw(4) << heterogen->GetSequenceNumber();
        else
            stream << right << setw(4) << " ";
        stream << right << setw(1) << heterogen->GetInsertionCode()
               << left << setw(2) << " ";
        if(heterogen->GetNumberOfHeterogenAtoms() != iNotSet)
            stream << right << setw(5) << heterogen->GetNumberOfHeterogenAtoms();
        else
            stream << right << setw(5) << " ";
        stream << left << setw(5) << " "
               << left << setw(40) << heterogen->GetDscr()
               << left << setw(10) << " "
               << endl;
    }
}

void PdbFile::ResolveHeterogenNameCard(std::ofstream& stream)
{
    PdbHeterogenNameCard::HeterogenNameMap heterogen_name_map = heterogens_name_->GetHeterogenNames();
    for(PdbHeterogenNameCard::HeterogenNameMap::iterator it = heterogen_name_map.begin(); it != heterogen_name_map.end(); it++)
    {
        PdbHeterogenName* heterogen_name = (*it).second;
        const int MAX_NAME_LENGTH_IN_LINE = 55;
        if((int)heterogen_name->GetHeterogenName().length() > MAX_NAME_LENGTH_IN_LINE)
        {
            stream << left << setw(6) << heterogens_name_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << right << setw(3) << heterogen_name->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << left << setw(55) << heterogen_name->GetHeterogenName().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << left << setw(10) << " "
                   << endl;
            int counter = ceil((double)(heterogen_name->GetHeterogenName().length()) / MAX_NAME_LENGTH_IN_LINE);
            for(int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << left << setw(6) << heterogens_name_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << right << setw(3) << heterogen_name->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << left << setw(55) << heterogen_name->GetHeterogenName().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),MAX_NAME_LENGTH_IN_LINE)
                           << left << setw(10) << " "
                           << endl;
                }
                else
                {
                    stream << left << setw(6) << heterogens_name_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << right << setw(3) << heterogen_name->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << left << setw(55) << heterogen_name->GetHeterogenName().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),heterogen_name->GetHeterogenName().length()-MAX_NAME_LENGTH_IN_LINE*(i-1))
                           << left << setw(10) << " "
                           << endl;
                }
            }
        }
        else
        {
            stream << left << setw(6) << heterogens_name_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << right << setw(3) << heterogen_name->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << left << setw(55) << heterogen_name->GetHeterogenName()
                   << left << setw(10) << " "
                   << endl;
        }
    }
}

void PdbFile::ResolveHeterogenSynonymCard(std::ofstream& stream)
{
    PdbHeterogenSynonymCard::HeterogenSynonymMap heterogen_synonym_map = heterogen_synonyms_->GetHeterogensSynonyms();
    for(PdbHeterogenSynonymCard::HeterogenSynonymMap::iterator it = heterogen_synonym_map.begin(); it != heterogen_synonym_map.end(); it++)
    {
        PdbHeterogenSynonym* heterogen_synonym = (*it).second;
        stringstream ss;
        vector<string> synonyms = heterogen_synonym->GetHeterogenSynonyms();
        for(vector<string>::iterator it = synonyms.begin(); it != synonyms.end(); it++)
        {
            if(it != synonyms.end() - 1)
            {
                ss << (*it) << "; ";
            }
            else
            {
                ss << (*it) << ";";
            }
        }
        const int MAX_SYNONYM_LENGTH_IN_LINE = 55;
        if((int)ss.str().length() > MAX_SYNONYM_LENGTH_IN_LINE)
        {
            stream << left << setw(6) << heterogen_synonyms_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << right << setw(3) << heterogen_synonym->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << left << setw(55) << ss.str().substr(0,MAX_SYNONYM_LENGTH_IN_LINE)
                   << left << setw(10) << " "
                   << endl;
            int counter = ceil((double)(ss.str().length()) / MAX_SYNONYM_LENGTH_IN_LINE);
            for(int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << left << setw(6) << heterogen_synonyms_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << right << setw(3) << heterogen_synonym->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << left << setw(55) << ss.str().substr(MAX_SYNONYM_LENGTH_IN_LINE*(i-1),MAX_SYNONYM_LENGTH_IN_LINE)
                           << left << setw(10) << " "
                           << endl;
                }
                else
                {
                    stream << left << setw(6) << heterogen_synonyms_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << right << setw(3) << heterogen_synonym->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << left << setw(55) << ss.str().substr(MAX_SYNONYM_LENGTH_IN_LINE*(i-1),ss.str().length()-MAX_SYNONYM_LENGTH_IN_LINE*(i-1))
                           << left << setw(10) << " "
                           << endl;
                }
            }
        }
        else
        {
            stream << left << setw(6) << heterogen_synonyms_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << right << setw(3) << heterogen_synonym->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << left << setw(55) << ss.str()
                   << left << setw(10) << " "
                   << endl;
        }
    }
}

void PdbFile::ResolveFormulaCard(std::ofstream& stream)
{
    PdbFormulaCard::FormulaMap formula_map = formulas_->GetFormulas();
    for(PdbFormulaCard::FormulaMap::iterator it = formula_map.begin(); it != formula_map.end(); it++)
    {
        PdbFormula* formula = (*it).second;
        const int MAX_NAME_LENGTH_IN_LINE = 51;
        if((int)formula->GetChemicalFormula().length() > MAX_NAME_LENGTH_IN_LINE)
        {
            stream << left << setw(6) << formulas_->GetRecordName()
                   << left << setw(2) << " ";
            if(formula->GetComponentNumber() != iNotSet)
                stream << right << setw(2) << formula->GetComponentNumber();
            else
                stream << right << setw(2) << " ";
            stream << left << setw(2) << " "
                   << right << setw(3) << formula->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << left << setw(51) << formula->GetChemicalFormula().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << left << setw(10) << " "
                   << endl;
            int counter = ceil((double)(formula->GetChemicalFormula().length()) / MAX_NAME_LENGTH_IN_LINE);
            for(int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << left << setw(6) << formulas_->GetRecordName()
                           << left << setw(2) << " ";
                    if(formula->GetComponentNumber() != iNotSet)
                        stream << right << setw(2) << formula->GetComponentNumber();
                    else
                        stream << right << setw(2) << " ";
                    stream << left << setw(2) << " "
                           << right << setw(3) << formula->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << left << setw(51) << formula->GetChemicalFormula().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),MAX_NAME_LENGTH_IN_LINE)
                           << left << setw(10) << " "
                           << endl;
                }
                else
                {
                    stream << left << setw(6) << formulas_->GetRecordName()
                           << left << setw(2) << " ";
                    if(formula->GetComponentNumber()!= iNotSet)
                        stream << right << setw(2) << formula->GetComponentNumber();
                    else
                        stream << right << setw(2) << " ";
                    stream << left << setw(2) << " "
                           << right << setw(3) << formula->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << left << setw(51) << formula->GetChemicalFormula().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),formula->GetChemicalFormula().length()-MAX_NAME_LENGTH_IN_LINE*(i-1))
                           << left << setw(10) << " "
                           << endl;
                }
            }
        }
        else
        {
            stream << left << setw(6) << formulas_->GetRecordName()
                   << left << setw(2) << " ";
            if(formula->GetComponentNumber() != iNotSet)
                stream << right << setw(2) << formula->GetComponentNumber();
            else
                stream << right << setw(2) << " ";
            stream << left << setw(2) << " "
                   << right << setw(3) << formula->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << left << setw(51) << formula->GetChemicalFormula().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << left << setw(10) << " "
                   << endl;
        }
    }
}

void PdbFile::ResolveHelixCard(std::ofstream& stream)
{
    PdbHelixCard::HelixMap helix_map = helixes_->GetHelixes();
    int counter = helix_map.size();
    int serial_number = 1;
    while(serial_number <= counter)
    {
        for(PdbHelixCard::HelixMap::iterator it = helix_map.begin(); it != helix_map.end(); it++)
        {
            
            PdbHelix* helix = (*it).second;
            PdbHelix::HelixResidueVector helix_residues = helix->GetHelixResidues();
            if(helix->GetHelixSerialNumber() == serial_number)
            {
                stream << left << setw(6) << helixes_->GetRecordName()
                       << left << setw(1) << " ";
                if(helix->GetHelixSerialNumber() != iNotSet)
                    stream << right << setw(3) << helix->GetHelixSerialNumber();
                else
                    stream << right << setw(3) << " ";
                stream << left << setw(1) << " "
                       << right << setw(3) << helix->GetHelixId()
                       << left << setw(1) << " "
                       << right << setw(3) << helix_residues.at(0)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << helix_residues.at(0)->GetResidueChainId()
                       << left << setw(1) << " ";
                if(helix_residues.at(0)->GetResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << helix_residues.at(0)->GetResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << helix_residues.at(0)->GetResidueInsertionCode()
                       << left << setw(1) << " "
                       << right << setw(3) << helix_residues.at(1)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << helix_residues.at(1)->GetResidueChainId()
                       << left << setw(1) << " ";
                if(helix_residues.at(1)->GetResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << helix_residues.at(1)->GetResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << helix_residues.at(1)->GetResidueInsertionCode();
                if(helix->GetHelixClass() != UnknownHelix)
                    stream << right << setw(2) << helix->GetHelixClass();
                else
                    stream << right << setw(2) << " ";
                stream << left << setw(30) << helix->GetComment()
                       << left << setw(1) << " ";
                if(helix->GetHelixLength() != dNotSet)
                    stream << right << setw(5) << helix->GetHelixLength();
                else
                    stream << right << setw(5) << " ";
                stream << left << setw(4) << " "
                       << endl;
                break;
            }
        }
        serial_number++;
    }
}

void PdbFile::ResolveSheetCard(std::ofstream& stream)
{
    PdbSheetCard::SheetMap sheet_map = sheets_->GetSheets();
    for(PdbSheetCard::SheetMap::iterator it = sheet_map.begin(); it != sheet_map.end(); it++)
    {
        PdbSheet* sheet = (*it).second;
        PdbSheet::SheetStrandVector strands = sheet->GetStrands();
        int serial_number = 1;
        for(PdbSheet::SheetStrandVector::iterator it1 = strands.begin(); it1 != strands.end(); it1++)
        {
            PdbSheetStrand* strand = (*it1);
            PdbSheetStrand::SheetStrandResidueVector strand_residues = strand->GetStrandResidues();
            if(strand->GetSense() != 0)
            {
                stream << left << setw(6) << sheets_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << serial_number
                       << left << setw(1) << " "
                       << right << setw(3) << sheet->GetSheetId();
                if(sheet->GetNumberOfStrands() != iNotSet)
                    stream << right << setw(2) << sheet->GetNumberOfStrands();
                else
                    stream << right << setw(2) << " ";
                stream << left << setw(1) << " "
                       << right << setw(3) << strand_residues.at(0)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(0)->GetResidueChainId();
                if(strand_residues.at(0)->GetResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << strand_residues.at(0)->GetResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << strand_residues.at(0)->GetResidueInsertionCode()
                       << left << setw(1) << " "
                       << right << setw(3) << strand_residues.at(1)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(1)->GetResidueChainId();
                if(strand_residues.at(1)->GetResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << strand_residues.at(1)->GetResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << strand_residues.at(1)->GetResidueInsertionCode();
                if(strand->GetSense() != UnknownStrand)
                    stream << right << setw(2) << strand->GetSense();
                else
                    stream << right << setw(2) << " ";
                stream << left << setw(1) << " "
                       << left << setw(4) << strand->GetCurrentAtom()
                       << right << setw(3) << strand_residues.at(2)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(2)->GetResidueChainId();
                if(strand_residues.at(2)->GetResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << strand_residues.at(2)->GetResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << strand_residues.at(2)->GetResidueInsertionCode()
                       << left << setw(1) << " "
                       << left << setw(4) << strand->GetPreviousAtom()
                       << right << setw(3) << strand_residues.at(3)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(3)->GetResidueChainId();
                if(strand_residues.at(3)->GetResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << strand_residues.at(3)->GetResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << strand_residues.at(3)->GetResidueInsertionCode()
                       << left << setw(10) << " "
                       << endl;
            }
            else
            {
                stream << left << setw(6) << sheets_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << serial_number
                       << left << setw(1) << " "
                       << right << setw(3) << sheet->GetSheetId();
                if(sheet->GetNumberOfStrands() != iNotSet)
                    stream << right << setw(2) << sheet->GetNumberOfStrands();
                else
                    stream << right << setw(2) << " ";
                stream << left << setw(1) << " "
                       << right << setw(3) << strand_residues.at(0)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(0)->GetResidueChainId();
                if(strand_residues.at(0)->GetResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << strand_residues.at(0)->GetResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << strand_residues.at(0)->GetResidueInsertionCode()
                       << left << setw(1) << " "
                       << right << setw(3) << strand_residues.at(1)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(1)->GetResidueChainId();
                if(strand_residues.at(1)->GetResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << strand_residues.at(1)->GetResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                stream << right << setw(1) << strand_residues.at(1)->GetResidueInsertionCode();
                if(strand->GetSense() != UnknownStrand)
                    stream << right << setw(2) << strand->GetSense();
                else
                    stream << right << setw(2) << " ";
                stream << left << setw(40) << " "
                       << endl;
            }
            serial_number++;
        }
        
    }
}

void PdbFile::ResolveDisulfideBondCard(std::ofstream& stream)
{
    PdbDisulfideBondCard::DisulfideResidueBondMap disulfide_bond_map = disulfide_bonds_->GetDisulfideResidueBonds();
    for(PdbDisulfideBondCard::DisulfideResidueBondMap::iterator it = disulfide_bond_map.begin(); it != disulfide_bond_map.end(); it++)
    {
        PdbDisulfideResidueBond* disulfide_bonds = (*it).second;
        PdbDisulfideResidueBond::DisulfideResidueVector disulfide_bonds_residues = disulfide_bonds->GetResidues();
        stream << left << setw(6) << disulfide_bonds_->GetRecordName()
               << left << setw(1) << " ";
        if(disulfide_bonds->GetSerialNumber() != iNotSet)
            stream << right << setw(3) << disulfide_bonds->GetSerialNumber();
        else
            stream << right << setw(3) << " ";
        stream << left << setw(1) << " "
               << right << setw(3) << disulfide_bonds_residues.at(0)->GetResidueName()
               << left << setw(1) << " "
               << right << setw(1) << disulfide_bonds_residues.at(0)->GetResidueChainId()
               << left << setw(1) << " ";
        if(disulfide_bonds_residues.at(0)->GetResidueSequenceNumber() != iNotSet)
            stream << right << setw(4) << disulfide_bonds_residues.at(0)->GetResidueSequenceNumber();
        else
            stream << right << setw(4) << " ";
        stream << right << setw(1) << disulfide_bonds_residues.at(0)->GetResidueInsertionCode()
               << left << setw(3) << " "
               << right << setw(3) << disulfide_bonds_residues.at(1)->GetResidueName()
               << left << setw(1) << " "
               << right << setw(1) << disulfide_bonds_residues.at(1)->GetResidueChainId()
               << left << setw(1) << " ";
        if(disulfide_bonds_residues.at(1)->GetResidueSequenceNumber() != iNotSet)
            stream << right << setw(4) << disulfide_bonds_residues.at(1)->GetResidueSequenceNumber();
        else
            stream << right << setw(4) << " ";
        stream << right << setw(1) << disulfide_bonds_residues.at(1)->GetResidueInsertionCode()
               << left << setw(23) << " ";
        if(disulfide_bonds_residues.at(0)->GetSymmetryOperator() != iNotSet)
            stream << right << setw(6) << disulfide_bonds_residues.at(0)->GetSymmetryOperator();
        else
            stream << right << setw(6) << " ";
        stream << left << setw(1) << " ";
        if(disulfide_bonds_residues.at(1)->GetSymmetryOperator() != iNotSet)
            stream << right << setw(6) << disulfide_bonds_residues.at(1)->GetSymmetryOperator();
        else
            stream << right << setw(6) << " ";
        stream << left << setw(1) << " ";
        if(disulfide_bonds->GetBondLength() != dNotSet)
            stream << right << setw(5) << fixed << setprecision(2) << disulfide_bonds->GetBondLength();
        else
            stream << right << setw(5) << " ";
        stream << left << setw(2) << " "
               << endl;
    }
}

void PdbFile::ResolveLinkCard(std::ofstream& stream)
{
    PdbLinkCard::LinkVector links = links_->GetResidueLinks();
    for(PdbLinkCard::LinkVector::iterator it = links.begin(); it != links.end(); it++)
    {
        PdbLink* link = (*it);
        PdbLink::LinkResidueVector link_residues = link->GetResidues();
        stream << left << setw(6) << links_->GetRecordName()
               << left << setw(6) << " "
               << left << setw(4) << link_residues.at(0)->GetAtomName();
        if(link_residues.at(0)->GetAlternateLocationIndicator() != BLANK_SPACE)
            stream << right << setw(1) << link_residues.at(0)->GetAlternateLocationIndicator();
        else
            stream << right << setw(1) << " ";
        stream << right << setw(3) << link_residues.at(0)->GetResidueName()
               << left << setw(1) << " ";
        if(link_residues.at(0)->GetResidueChainId() != BLANK_SPACE)
            stream << right << setw(1) << link_residues.at(0)->GetResidueChainId();
        else
            stream << right << setw(1) << " ";
        if(link_residues.at(0)->GetResidueSequenceNumber() != iNotSet)
            stream << right << setw(4) << link_residues.at(0)->GetResidueSequenceNumber();
        else
            stream << right << setw(4) << " ";
        if(link_residues.at(0)->GetResidueInsertionCode() != BLANK_SPACE)
            stream << right << setw(1) << link_residues.at(0)->GetResidueInsertionCode();
        else
            stream << right << setw(1) << " ";

        stream << left << setw(15) << " "
               << left << setw(4) << link_residues.at(1)->GetAtomName();
        if(link_residues.at(1)->GetAlternateLocationIndicator() != BLANK_SPACE)
            stream << right << setw(1) << link_residues.at(1)->GetAlternateLocationIndicator();
        else
            stream << right << setw(1) << " ";
        stream << right << setw(3) << link_residues.at(1)->GetResidueName()
               << left << setw(1) << " ";
        if(link_residues.at(1)->GetResidueChainId() != BLANK_SPACE)
            stream << right << setw(1) << link_residues.at(1)->GetResidueChainId();
        else
            stream << right << setw(1) << " ";
        if(link_residues.at(1)->GetResidueSequenceNumber() != iNotSet)
            stream << right << setw(4) << link_residues.at(1)->GetResidueSequenceNumber();
        else
            stream << right << setw(4) << " ";
        if(link_residues.at(1)->GetResidueInsertionCode() != BLANK_SPACE)
            stream << right << setw(1) << link_residues.at(1)->GetResidueInsertionCode();
        else
            stream << right << setw(1) << " ";
        stream << left << setw(2) << " ";
        if(link_residues.at(0)->GetSymmetryOperator() != iNotSet)
            stream << right << setw(6) << link_residues.at(0)->GetSymmetryOperator();
        else
            stream << right << setw(6) << " ";
        stream << left << setw(1) << " ";
        if(link_residues.at(1)->GetSymmetryOperator() != iNotSet)
            stream << right << setw(6) << link_residues.at(1)->GetSymmetryOperator();
        else
            stream << right << setw(6) << " ";
        stream << left << setw(1) << " ";
        if(link->GetLinkLength() != dNotSet)
            stream << right << setw(5) << fixed << setprecision(2) << link->GetLinkLength();
        else
            stream << right << setw(5) << " ";
        stream << left << setw(2) << " "
               << endl;
    }
}

void PdbFile::ResolveCISPeptideCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveSiteCard(std::ofstream& stream)
{
    PdbSiteCard::PdbSiteMap site_map = sites_->GetResidueSites();
    for(PdbSiteCard::PdbSiteMap::iterator it = site_map.begin(); it != site_map.end(); it++)
    {
        PdbSite* site = (*it).second;
        PdbSite::SiteResidueVector site_residues = site->GetResidues();
        const int MAX_RESIDUE_IN_LINE = 4;
        const int RESIDUE_LENGHT_IN_LINE = 11;
        int number_of_residues = site_residues.size();
        
        int sequence_number = 1;
        if(number_of_residues > MAX_RESIDUE_IN_LINE)
        {
            int number_of_lines = ceil((double)(number_of_residues) / MAX_RESIDUE_IN_LINE);
            while(sequence_number <= number_of_lines)
            {
                if(sequence_number != number_of_lines)
                {
                    stringstream ss;
                    for(PdbSite::SiteResidueVector::iterator it1 = site_residues.begin()+(sequence_number - 1)*MAX_RESIDUE_IN_LINE;
                        it1 != site_residues.begin()+(sequence_number)*MAX_RESIDUE_IN_LINE; it1++)
                    {
                        PdbSiteResidue* residue = (*it1);
                        ss << left << setw(1) << " "
                           << right << setw(3) << residue->GetResidueName()
                           << left << setw(1) << " "
                           << right << setw(1) << residue->GetResidueChainId();
                        if(residue->GetResidueSequenceNumber() != iNotSet)
                            ss << right << setw(4) << residue->GetResidueSequenceNumber();
                        else
                            ss << right << setw(4) << " ";
                        ss << right << setw(1) << residue->GetResidueInsertionCode();
                    }
                    ss << left << setw(19) << " ";
                    stream << left << setw(6) << sites_->GetRecordName()
                           << left << setw(1) << " "
                           << right << setw(3) << sequence_number
                           << left << setw(1) << " "
                           << right << setw(3) << site->GetSiteId()
                           << left << setw(1) << " ";
                    if(site->GetNumberOfResidues() != iNotSet)
                        stream << right << setw(2) << site->GetNumberOfResidues();
                    else
                        stream << right << setw(2) << " ";
                    stream << left << setw(63) << ss.str()
                           << endl;
                }
                else
                {
                    stringstream ss;
                    for(PdbSite::SiteResidueVector::iterator it1 = site_residues.begin()+(sequence_number - 1)*MAX_RESIDUE_IN_LINE;
                        it1 != site_residues.end(); it1++)
                    {
                        PdbSiteResidue* residue = (*it1);
                        ss << left << setw(1) << " "
                           << right << setw(3) << residue->GetResidueName()
                           << left << setw(1) << " "
                           << right << setw(1) << residue->GetResidueChainId();
                        if(residue->GetResidueSequenceNumber() != iNotSet)
                            ss << right << setw(4) << residue->GetResidueSequenceNumber();
                        else
                            ss << right << setw(4) << " ";
                        ss << right << setw(1) << residue->GetResidueInsertionCode();
                    }
                    if((sequence_number*MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE != 0)
                        ss << left << setw((sequence_number*MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE) << " ";
                    ss << left << setw(19) << " ";
                    stream << left << setw(6) << sites_->GetRecordName()
                           << left << setw(1) << " "
                           << right << setw(3) << sequence_number
                           << left << setw(1) << " "
                           << right << setw(3) << site->GetSiteId()
                           << left << setw(1) << " ";
                    if(site->GetNumberOfResidues() != iNotSet)
                        stream << right << setw(2) << site->GetNumberOfResidues();
                    else
                        stream << right << setw(2) << " ";
                    stream << left << setw(63) << ss.str()
                           << endl;
                }
                sequence_number++;
            }
            
        }
        else
        {
            stringstream ss;
            for(PdbSite::SiteResidueVector::iterator it1 = site_residues.begin(); it1 != site_residues.end(); it1++)
            {
                PdbSiteResidue* residue = (*it1);
                ss << left << setw(1) << " "
                   << right << setw(3) << residue->GetResidueName()
                   << left << setw(1) << " "
                   << right << setw(1) << residue->GetResidueChainId();
                if(residue->GetResidueSequenceNumber() != iNotSet)
                    ss << right << setw(4) << residue->GetResidueSequenceNumber();
                else
                    ss << right << setw(4) << " ";
                ss << right << setw(1) << residue->GetResidueInsertionCode();
            }
            if((MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE != 0)
                ss << left << setw((MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE) << " ";
            ss << left << setw(19) << " ";
            
            stream << left << setw(6) << sites_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << sequence_number
                   << left << setw(1) << " "
                   << right << setw(3) << site->GetSiteId()
                   << left << setw(1) << " ";
            if(site->GetNumberOfResidues() != iNotSet)
                stream << right << setw(2) << site->GetNumberOfResidues();
            else
                stream << right << setw(2) << " ";
            stream << left << setw(63) << ss.str()
                   << endl;
            sequence_number++;
        }
    }
}

void PdbFile::ResolveCrystallographyCard(std::ofstream& stream)
{
    stream << left << setw(6) << crystallography_->GetRecordName();
    if(crystallography_->GetA() != dNotSet)
        stream << right << setw(9) << fixed << setprecision(3) << crystallography_->GetA();
    else
        stream << right << setw(9) << " ";
    if(crystallography_->GetB() != dNotSet)
        stream << right << setw(9) << fixed << setprecision(3) << crystallography_->GetB();
    else
        stream << right << setw(9) << " ";
    if(crystallography_->GetC() != dNotSet)
        stream << right << setw(9) << fixed << setprecision(3) << crystallography_->GetC();
    else
        stream << right << setw(9) << " ";
    if(crystallography_->GetAlpha() != dNotSet)
        stream << right << setw(7) << fixed << setprecision(2) << crystallography_->GetAlpha();
    else
        stream << right << setw(7) << " ";
    if(crystallography_->GetBeta() != dNotSet)
        stream << right << setw(7) << fixed << setprecision(2) << crystallography_->GetBeta();
    else
        stream << right << setw(7) << " ";
    if(crystallography_->GetGamma() != dNotSet)
        stream << right << setw(7) << fixed << setprecision(2) << crystallography_->GetGamma();
    else
        stream << right << setw(7) << " ";
    stream << left << setw(1) << " "
           << left << setw(11) << crystallography_->GetSpaceGroup();
    if(crystallography_->GetZValue() != iNotSet)
        stream << right << setw(4) << crystallography_->GetZValue();
    else
        stream << right << setw(4) << " ";
    stream << left << setw(10) << " "
           << endl;
}

void PdbFile::ResolveOriginCard(std::ofstream& stream)
{
    PdbOriginXnCard::OriginXnVector origins = origins_->GetOriginXN();
    for(PdbOriginXnCard::OriginXnVector::iterator it = origins.begin(); it != origins.end(); it++)
    {
        PdbOriginXn* origin = (*it);
        stringstream ss;
        ss << origin->GetRecordName() << origin->GetN();
        stream << left << setw(6) << ss.str()
               << left << setw(4) << " ";
        if(origin->GetOrigin().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
        {
            stream << right << setw(10) << fixed << setprecision(6) << origin->GetOrigin().GetX()
                   << right << setw(10) << fixed << setprecision(6) << origin->GetOrigin().GetY()
                   << right << setw(10) << fixed << setprecision(6) << origin->GetOrigin().GetZ();
        }
        else
        {
            stream << right << setw(10) << " "
                   << right << setw(10) << " "
                   << right << setw(10) << " ";
        }
        
        stream << left << setw(5) << " ";
        if(origin->GetT() != dNotSet)
            stream << right << setw(10) << fixed << setprecision(5) << origin->GetT();
        else
            stream << right << setw(10) << " ";
        stream << left << setw(25) << " "
               << endl;
    }
}

void PdbFile::ResolveScaleCard(std::ofstream& stream)
{
    PdbScaleNCard::ScaleNVector scales = scales_->GetScaleN();
    for(PdbScaleNCard::ScaleNVector::iterator it = scales.begin(); it != scales.end(); it++)
    {
        PdbScaleN* scale = (*it);
        stringstream ss;
        ss << scale->GetRecordName() << scale->GetN();
        stream << left << setw(6) << ss.str()
               << left << setw(4) << " ";
        if(scale->GetScaleVector().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
        {
            stream << right << setw(10) << fixed << setprecision(6) << scale->GetScaleVector().GetX()
                   << right << setw(10) << fixed << setprecision(6) << scale->GetScaleVector().GetY()
                   << right << setw(10) << fixed << setprecision(6) << scale->GetScaleVector().GetZ();
        }
        else
        {
            stream << right << setw(10) << " "
                   << right << setw(10) << " "
                   << right << setw(10) << " ";
        }
        stream << left << setw(5) << " ";
        if(scale->GetU() != dNotSet)
            stream << right << setw(10) << fixed << setprecision(5) << scale->GetU();
        else
            stream << right << setw(10) << " ";
        stream << left << setw(25) << " "
               << endl;
    }
}

void PdbFile::ResolveMatrixCard(std::ofstream& stream)
{
    PdbMatrixNCard::MatrixNVectorVector matrices = matrices_->GetMatrixN();
    int number_of_matrix_entries = matrices.at(0).size();
    for(int i = 0; i < number_of_matrix_entries; i++)
    {
        for(unsigned int j = 0; j < 3; j++)
        {
            PdbMatrixNCard::MatrixNVector matrix_vector = matrices.at(j);
            PdbMatrixN* matrix = matrix_vector.at(i);
            stringstream ss;
            ss << matrix->GetRecordName() << matrix->GetN();
            stream << left << setw(6) << ss.str()
                   << left << setw(1) << " ";
            if(matrix->GetSerialNumber() != iNotSet)
                stream << right << setw(3) << matrix->GetSerialNumber();
            else
                stream << right << setw(3) << " ";
            if(matrix->GetTransformationVector().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
            {
                stream << right << setw(10) << fixed << setprecision(6) << matrix->GetTransformationVector().GetX()
                       << right << setw(10) << fixed << setprecision(6) << matrix->GetTransformationVector().GetY()
                       << right << setw(10) << fixed << setprecision(6) << matrix->GetTransformationVector().GetZ();
            }
            else
            {
                stream << right << setw(10) << " "
                       << right << setw(10) << " "
                       << right << setw(10) << " ";
            }
            stream << left << setw(5) << " ";
            if(matrix->GetV() != dNotSet)
                stream << right << setw(10) << fixed << setprecision(5) << matrix->GetV();
            else
                stream << right << setw(10) << " ";
            stream << left << setw(4) << " ";
            if(matrix->GetIGiven() != iNotSet)
                stream << right << setw(1) << matrix->GetIGiven();
            else
                stream << right << setw(1) << " ";
            stream << left << setw(20) << " "
                   << endl;
        }
    }
    
}

void PdbFile::ResolveModelCard(std::ofstream& stream)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    int number_of_models = models.size();
    if(number_of_models == 1)
    {
        for(PdbModelCard::PdbModelMap::iterator it = models.begin(); it != models.end(); it++)
        {
            PdbModel* model = (*it).second;
            PdbModelResidueSet* residue_set = model->GetModelResidueSet();
            PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
            for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
            {
                PdbAtomCard* atom_card = (*it1);
                int serial_number = 0;
                string residue_name = "";
                char chain_id = ' ';
                int residue_sequence_number = 0;
                char insertion_code = ' ';
                PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
                int atoms_size = ordered_atoms.size();
                for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
                {
                    PdbAtom* atom = (*it2);
                    stream << left << setw(6) << atom_card->GetRecordName();
                    if(atom->GetAtomSerialNumber() != iNotSet)
                        stream << right << setw(5) << atom->GetAtomSerialNumber();
                    else
                        stream << right << setw(5) << " ";
                    stream << left << setw(1) << " "
                           << left << setw(4) << atom->GetAtomName();
                    if(atom->GetAtomAlternateLocation() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << atom->GetAtomAlternateLocation();
                    stream << right << setw(3) << atom->GetAtomResidueName()
                           << left << setw(1) << " ";
                    if(atom->GetAtomChainId() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << atom->GetAtomChainId();
                    if(atom->GetAtomResidueSequenceNumber() != iNotSet)
                        stream << right << setw(4) << atom->GetAtomResidueSequenceNumber();
                    else
                        stream << right << setw(4) << " ";
                    if(atom->GetAtomInsertionCode() == BLANK_SPACE)
                        stream << left << setw(1) <<  ' ';
                    else
                        stream << left << setw(1) << atom->GetAtomInsertionCode();
                    stream << left << setw(3) << " ";
                    if(atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
                    {
                        stream << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetX()
                               << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetY()
                               << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetZ();
                    }
                    else
                    {
                        stream << right << setw(8) << " "
                               << right << setw(8) << " "
                               << right << setw(8) << " ";
                    }
                    if(atom->GetAtomOccupancy() != dNotSet)
                        stream << right << setw(6) << fixed << setprecision(2) << atom->GetAtomOccupancy();
                    else
                        stream << right << setw(6) << " ";
                    if(atom->GetAtomTempretureFactor() != dNotSet)
                        stream << right << setw(6) << fixed << setprecision(2) << atom->GetAtomTempretureFactor();
                    else
                        stream << right << setw(6) << " ";
                    stream << left << setw(10) << " "
                           << right << setw(2) << atom->GetAtomElementSymbol()
                           << left << setw(2) << atom->GetAtomCharge()
                           << endl;
                    serial_number = atom->GetAtomSerialNumber();
                    residue_name = atom->GetAtomResidueName();
                    chain_id = atom->GetAtomChainId();
                    residue_sequence_number = atom->GetAtomResidueSequenceNumber();
                }
                if(atoms_size != 0)
                {
                    stream << left << setw(6) << "TER";
                    if(serial_number != iNotSet)
                        stream << right << setw(5) << (serial_number+1);
                    else
                        stream << right << setw(5) << " ";
                    stream << left << setw(6) << " "
                           << right << setw(3) << residue_name
                           << left << setw(1) << " ";
                    if(chain_id == BLANK_SPACE)
                        stream << left << setw(1) << " ";
                    else
                        stream << left << setw(1) << chain_id;
                    if(residue_sequence_number != iNotSet)
                        stream << right << setw(4) << residue_sequence_number;
                    else
                        stream << right << setw(4) << " ";
                    stream << left << setw(1) << insertion_code
                           << left << setw(53) << " "
                           << endl;
                }
            }
            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
            for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
            {
                PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
                PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
                for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
                {
                    PdbAtom* heterogen_atom = (*it2);
                    stream << left << setw(6) << heterogen_atom_card->GetRecordName();
                    if(heterogen_atom->GetAtomSerialNumber() != iNotSet)
                        stream << right << setw(5) << heterogen_atom->GetAtomSerialNumber();
                    else
                        stream << right << setw(5) << " ";
                    stream << left << setw(1) << " "
                           << left << setw(4) << heterogen_atom->GetAtomName();
                    if(heterogen_atom->GetAtomAlternateLocation() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << heterogen_atom->GetAtomAlternateLocation();
                    stream << right << setw(3) << heterogen_atom->GetAtomResidueName()
                           << left << setw(1) << " ";
                    if(heterogen_atom->GetAtomChainId() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << heterogen_atom->GetAtomChainId();
                    if(heterogen_atom->GetAtomResidueSequenceNumber() != iNotSet)
                        stream << right << setw(4) << heterogen_atom->GetAtomResidueSequenceNumber();
                    else
                        stream << right << setw(4) << " ";
                    if(heterogen_atom->GetAtomInsertionCode() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << heterogen_atom->GetAtomInsertionCode();
                    stream << left << setw(3) << " ";
                    if(heterogen_atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
                    {
                        stream << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetX()
                               << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetY()
                               << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetZ();
                    }
                    else
                    {
                        stream << right << setw(8) << " "
                               << right << setw(8) << " "
                               << right << setw(8) << " ";
                    }
                    if(heterogen_atom->GetAtomOccupancy() != dNotSet)
                        stream << right << setw(6) << fixed << setprecision(2) << heterogen_atom->GetAtomOccupancy();
                    else
                        stream << right << setw(6) << " ";
                    if(heterogen_atom->GetAtomTempretureFactor() != dNotSet)
                        stream << right << setw(6) << fixed << setprecision(2) << heterogen_atom->GetAtomTempretureFactor();
                    else
                        stream << right << setw(6) << " ";
                    stream << left << setw(10) << " "
                           << right << setw(2) << heterogen_atom->GetAtomElementSymbol()
                           << left << setw(2) << heterogen_atom->GetAtomCharge()
                           << endl;
                }
            }
        }
    }
    else
    {
        for(PdbModelCard::PdbModelMap::iterator it = models.begin(); it != models.end(); it++)
        {
            PdbModel* model = (*it).second;
            stream << left << setw(6) << models_->GetRecordName()
                   << left << setw(4) << " ";
            if(model->GetModelSerialNumber() != iNotSet)
                stream << right << setw(4) << model->GetModelSerialNumber();
            else
                stream << right << setw(4) << " ";
            stream << left << setw(66) << " "
                   << endl;
            PdbModelResidueSet* residue_set = model->GetModelResidueSet();
            PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
            for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
            {
                PdbAtomCard* atom_card = (*it1);
                int serial_number = 0;
                string residue_name = "";
                char chain_id = ' ';
                int residue_sequence_number = 0;
                char insertion_code = ' ';
                PdbAtomCard::PdbAtomOrderVector ordered_atoms = atom_card->GetOrderedAtoms();
                int atoms_size = ordered_atoms.size();
                for(PdbAtomCard::PdbAtomOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
                {
                    PdbAtom* atom = (*it2);
                    stream << left << setw(6) << atom_card->GetRecordName();
                    if(atom->GetAtomSerialNumber() != iNotSet)
                        stream << right << setw(5) << atom->GetAtomSerialNumber();
                    else
                        stream << right << setw(5) << " ";
                    stream << left << setw(1) << " "
                           << left << setw(4) << atom->GetAtomName();
                    if(atom->GetAtomAlternateLocation() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << atom->GetAtomAlternateLocation();
                    stream << right << setw(3) << atom->GetAtomResidueName()
                           << left << setw(1) << " ";
                    if(atom->GetAtomChainId() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << atom->GetAtomChainId();
                    if(atom->GetAtomResidueSequenceNumber() != iNotSet)
                        stream << right << setw(4) << atom->GetAtomResidueSequenceNumber();
                    else
                        stream << right << setw(4) << " ";
                    if(atom->GetAtomInsertionCode() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << atom->GetAtomInsertionCode();
                    stream << left << setw(3) << " ";
                    if(atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
                        stream << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetX()
                               << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetY()
                               << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetZ();
                    else
                        stream << right << setw(8) << " "
                               << right << setw(8) << " "
                               << right << setw(8) << " ";
                    if(atom->GetAtomOccupancy() != dNotSet)
                        stream << right << setw(6) << fixed << setprecision(2) << atom->GetAtomOccupancy();
                    else
                        stream << right << setw(6) << " ";
                    if(atom->GetAtomTempretureFactor() != dNotSet)
                        stream << right << setw(6) << fixed << setprecision(2) << atom->GetAtomTempretureFactor();
                    else
                        stream << right << setw(6) << " ";
                    stream << left << setw(10) << " "
                           << right << setw(2) << atom->GetAtomElementSymbol()
                           << left << setw(2) << atom->GetAtomCharge()
                           << endl;
                    serial_number = atom->GetAtomSerialNumber();
                    residue_name = atom->GetAtomResidueName();
                    chain_id = atom->GetAtomChainId();
                    residue_sequence_number = atom->GetAtomResidueSequenceNumber();
                }
                if(atoms_size != 0)
                {
                    stream << left << setw(6) << "TER";
                    if(serial_number != iNotSet)
                        stream << right << setw(5) << (serial_number+1);
                    else
                        stream << right << setw(5) << " ";
                    stream << left << setw(6) << " "
                           << right << setw(3) << residue_name
                           << left << setw(1) << " ";
                    if(chain_id == BLANK_SPACE)
                        stream << left << setw(1) << " ";
                    else
                        stream << left << setw(1) << chain_id;
                    if(residue_sequence_number != iNotSet)
                        stream << right << setw(4) << residue_sequence_number;
                    else
                        stream << right << setw(4) << " ";
                    stream << left << setw(1) << insertion_code
                           << left << setw(53) << " "
                           << endl;
                }
            }
            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
            for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
            {
                PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
                PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtoms();
                for(PdbHeterogenAtomCard::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
                {
                    PdbAtom* heterogen_atom = (*it2);
                    stream << left << setw(6) << heterogen_atom_card->GetRecordName();
                    if(heterogen_atom->GetAtomSerialNumber() != iNotSet)
                        stream << right << setw(5) << heterogen_atom->GetAtomSerialNumber();
                    else
                        stream << right << setw(5) << " ";
                    stream << left << setw(1) << " "
                           << left << setw(4) << heterogen_atom->GetAtomName();
                    if(heterogen_atom->GetAtomAlternateLocation() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << heterogen_atom->GetAtomAlternateLocation();
                    stream << right << setw(3) << heterogen_atom->GetAtomResidueName()
                           << left << setw(1) << " ";
                    if(heterogen_atom->GetAtomChainId() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << heterogen_atom->GetAtomChainId();
                    if(heterogen_atom->GetAtomResidueSequenceNumber() != iNotSet)
                        stream << right << setw(4) << heterogen_atom->GetAtomResidueSequenceNumber();
                    else
                        stream << right << setw(4) << " ";
                    if(heterogen_atom->GetAtomInsertionCode() == BLANK_SPACE)
                        stream << left << setw(1) << ' ';
                    else
                        stream << left << setw(1) << heterogen_atom->GetAtomInsertionCode();
                    stream << left << setw(3) << " ";
                    if(heterogen_atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
                        stream << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetX()
                               << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetY()
                               << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetZ();
                    else
                        stream << right << setw(8) << " "
                               << right << setw(8) << " "
                               << right << setw(8) << " ";
                    if(heterogen_atom->GetAtomOccupancy() != dNotSet)
                        stream << right << setw(6) << fixed << setprecision(2) << heterogen_atom->GetAtomOccupancy();
                    else
                        stream << right << setw(6) << " ";
                    if(heterogen_atom->GetAtomTempretureFactor() != dNotSet)
                        stream << right << setw(6) << fixed << setprecision(2) << heterogen_atom->GetAtomTempretureFactor();
                    else
                        stream << right << setw(6) << " ";
                    stream << left << setw(10) << " "
                           << right << setw(2) << heterogen_atom->GetAtomElementSymbol()
                           << left << setw(2) << heterogen_atom->GetAtomCharge()
                           << endl;
                }
            }
            stream << left << setw(6) << "ENDMDL"
                   << left << setw(74) << " "
                   << endl;
        }
    }
}

void PdbFile::ResolveModelCardWithTheGivenModelNumber(std::ofstream& stream, int model_number)
{
    PdbModelCard::PdbModelMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModel* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtoms();
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomCard* atom_card = (*it1);
            int serial_number = 0;
            string residue_name = "";
            char chain_id = ' ';
            int residue_sequence_number = 0;
            char insertion_code = ' ';
            PdbAtomCard::PdbAtomMap atoms = atom_card->GetAtoms();
            int atoms_size = atoms.size();
            for(PdbAtomCard::PdbAtomMap::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                PdbAtom* atom = (*it2).second;
                stream << left << setw(6) << atom_card->GetRecordName();
                if(atom->GetAtomSerialNumber() != iNotSet)
                    stream << right << setw(5) << atom->GetAtomSerialNumber();
                else
                    stream << right << setw(5) << " ";
                stream << left << setw(1) << " "
                       << left << setw(4) << atom->GetAtomName();
                if(atom->GetAtomAlternateLocation() == BLANK_SPACE)
                    stream << left << setw(1) << ' ';
                else
                    stream << left << setw(1) << atom->GetAtomAlternateLocation();
                stream << right << setw(3) << atom->GetAtomResidueName()
                       << left << setw(1) << " ";
                if(atom->GetAtomChainId() == BLANK_SPACE)
                    stream << left << setw(1) << ' ';
                else
                    stream << left << setw(1) << atom->GetAtomChainId();
                if(atom->GetAtomResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << atom->GetAtomResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                if(atom->GetAtomInsertionCode() == BLANK_SPACE)
                    stream << left << setw(1) << ' ';
                else
                    stream << left << setw(1) << atom->GetAtomInsertionCode();
                stream << left << setw(3) << " ";
                if(atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
                {
                    stream << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetX()
                           << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetY()
                           << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetZ();
                }
                else
                {
                    stream << right << setw(8) << " "
                           << right << setw(8) << " "
                           << right << setw(8) << " ";
                }
                if(atom->GetAtomOccupancy() != dNotSet)
                    stream << right << setw(6) << fixed << setprecision(2) << atom->GetAtomOccupancy();
                else
                    stream << right << setw(6) << " ";
                if(atom->GetAtomTempretureFactor() != dNotSet)
                    stream << right << setw(6) << fixed << setprecision(2) << atom->GetAtomTempretureFactor();
                else
                    stream << right << setw(6) << " ";
                stream << left << setw(10) << " "
                       << right << setw(2) << atom->GetAtomElementSymbol()
                       << left << setw(2) << atom->GetAtomCharge()
                       << endl;
                serial_number = atom->GetAtomSerialNumber();
                residue_name = atom->GetAtomResidueName();
                chain_id = atom->GetAtomChainId();
                residue_sequence_number = atom->GetAtomResidueSequenceNumber();
            }
            if(atoms_size != 0)
            {
                stream << left << setw(6) << "TER";
                if(serial_number != iNotSet)
                    stream << right << setw(5) << (serial_number+1);
                else
                    stream << right << setw(5) << " ";
                stream << left << setw(6) << " "
                       << right << setw(3) << residue_name
                       << left << setw(1) << " ";
                if(chain_id == BLANK_SPACE)
                    stream << left << setw(1) << " ";
                else
                    stream << left << setw(1) << chain_id;
                if(residue_sequence_number != iNotSet)
                    stream << right << setw(4) << residue_sequence_number;
                else
                    stream << right << setw(4) << " ";
                stream << left << setw(1) << insertion_code
                       << left << setw(53) << " "
                       << endl;
            }
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtoms();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomCard* heterogen_atom_card = (*it1);
            PdbHeterogenAtomCard::PdbHeterogenAtomMap heterogen_atoms = heterogen_atom_card->GetHeterogenAtoms();
            for(PdbHeterogenAtomCard::PdbHeterogenAtomMap::iterator it2 = heterogen_atoms.begin(); it2 != heterogen_atoms.end(); it2++)
            {
                PdbAtom* heterogen_atom = (*it2).second;
                stream << left << setw(6) << heterogen_atom_card->GetRecordName();
                if(heterogen_atom->GetAtomSerialNumber() != iNotSet)
                    stream << right << setw(5) << heterogen_atom->GetAtomSerialNumber();
                else
                    stream << right << setw(5) << " ";
                stream << left << setw(1) << " "
                       << left << setw(4) << heterogen_atom->GetAtomName();
                if(heterogen_atom->GetAtomAlternateLocation() == BLANK_SPACE)
                    stream << left << setw(1) << ' ';
                else
                    stream << left << setw(1) << heterogen_atom->GetAtomAlternateLocation();
                stream << right << setw(3) << heterogen_atom->GetAtomResidueName()
                       << left << setw(1) << " ";
                if(heterogen_atom->GetAtomChainId() == BLANK_SPACE)
                    stream << left << setw(1) << ' ';
                else
                    stream << left << setw(1) << heterogen_atom->GetAtomChainId();
                if(heterogen_atom->GetAtomResidueSequenceNumber() != iNotSet)
                    stream << right << setw(4) << heterogen_atom->GetAtomResidueSequenceNumber();
                else
                    stream << right << setw(4) << " ";
                if(heterogen_atom->GetAtomInsertionCode() == BLANK_SPACE)
                    stream << left << setw(1) << ' ';
                else
                    stream << left << setw(1) << heterogen_atom->GetAtomInsertionCode();
                stream << left << setw(3) << " ";
                if(heterogen_atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
                {
                    stream << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetX()
                           << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetY()
                           << right << setw(8) << fixed << setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetZ();
                }
                else
                {
                    stream << right << setw(8) << " "
                           << right << setw(8) << " "
                           << right << setw(8) << " ";
                }
                if(heterogen_atom->GetAtomOccupancy() != dNotSet)
                    stream << right << setw(6) << fixed << setprecision(2) << heterogen_atom->GetAtomOccupancy();
                else
                    stream << right << setw(6) << " ";
                if(heterogen_atom->GetAtomTempretureFactor() != dNotSet)
                    stream << right << setw(6) << fixed << setprecision(2) << heterogen_atom->GetAtomTempretureFactor();
                else
                    stream << right << setw(6) << " ";
                stream << left << setw(10) << " "
                       << right << setw(2) << heterogen_atom->GetAtomElementSymbol()
                       << left << setw(2) << heterogen_atom->GetAtomCharge()
                       << endl;
            }
        }
    }
}

void PdbFile::ResolveConnectivityCard(std::ofstream& stream)
{
    PdbConnectCard::BondedAtomsSerialNumbersMap bonded_atoms = connectivities_->GetBondedAtomsSerialNumbers();
    for(PdbConnectCard::BondedAtomsSerialNumbersMap::iterator it = bonded_atoms.begin(); it != bonded_atoms.end(); it++)
    {
        vector<int> bonded_atoms_serial_number = (*it).second;
        int source_atom_serial_number = (*it).first;
        int number_of_bonded_atoms = bonded_atoms_serial_number.size();
        const int MAX_SERIAL_NUMBER_IN_LINE = 4;
        const int SERIAL_NUMBER_LENGTH = 5;
        if(number_of_bonded_atoms <= MAX_SERIAL_NUMBER_IN_LINE)
        {
            stream << left << setw(6) << connectivities_->GetRecordName();
            if(source_atom_serial_number != iNotSet)
                stream << right << setw(5) << source_atom_serial_number;
            else
                stream << right << setw(5) << " ";
            for(vector<int>::iterator it1 = bonded_atoms_serial_number.begin(); it1 != bonded_atoms_serial_number.end(); it1++)
            {
                int serial_number = (*it1);
                if(serial_number != iNotSet)
                    stream << right << setw(5) << serial_number;
                else
                    stream << right << setw(5) << " ";
            }
            if((MAX_SERIAL_NUMBER_IN_LINE-number_of_bonded_atoms)*SERIAL_NUMBER_LENGTH != 0)
                stream << left << setw((MAX_SERIAL_NUMBER_IN_LINE-number_of_bonded_atoms)*SERIAL_NUMBER_LENGTH) << " "
                       << left << setw(49) << " "
                       << endl;
            else
                stream << left << setw(49) << " "
                       << endl;
        }
        else
        {
            int number_of_lines = ceil((double)(number_of_bonded_atoms) / MAX_SERIAL_NUMBER_IN_LINE);
            for(int i = 1; i <= number_of_lines; i++)
            {
                stream << left << setw(6) << connectivities_->GetRecordName();
                if(source_atom_serial_number != iNotSet)
                    stream << right << setw(5) << source_atom_serial_number;
                else
                    stream << right << setw(5) << " ";
                if(i != number_of_lines)
                {
                    for(vector<int>::iterator it1 = bonded_atoms_serial_number.begin() + (i-1) * MAX_SERIAL_NUMBER_IN_LINE;
                        it1 != bonded_atoms_serial_number.begin() + i * MAX_SERIAL_NUMBER_IN_LINE; it1++)
                    {
                        int serial_number = (*it1);
                        if(serial_number != iNotSet)
                            stream << right << setw(5) << serial_number;
                        else
                            stream << right << setw(5) << " ";
                    }
                    stream << left << setw(49) << " "
                           << endl;
                }
                else
                {
                    for(vector<int>::iterator it1 = bonded_atoms_serial_number.begin() + (i-1) * MAX_SERIAL_NUMBER_IN_LINE;
                        it1 != bonded_atoms_serial_number.end(); it1++)
                    {
                        int serial_number = (*it1);
                        if(serial_number != iNotSet)
                            stream << right << setw(5) << serial_number;
                        else
                            stream << right << setw(5) << " ";
                    }
                    if((MAX_SERIAL_NUMBER_IN_LINE-(number_of_bonded_atoms-(i-1)*MAX_SERIAL_NUMBER_IN_LINE))*SERIAL_NUMBER_LENGTH != 0)
                        stream << left << setw((MAX_SERIAL_NUMBER_IN_LINE-number_of_bonded_atoms)*SERIAL_NUMBER_LENGTH) << " "
                               << left << setw(49) << " "
                               << endl;
                    else
                        stream << left << setw(49) << " "
                               << endl;
                }
            }
        }
    }
}

void PdbFile::ResolveMasterCard(std::ofstream& stream)
{
    
}

void PdbFile::ResolveEndCard(std::ofstream& stream)
{
    stream << left << setw(6) << "END" << left << setw(74) << " " << endl;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

void PdbFile::Print(ostream &out)
{
    if(header_ != NULL)
    {
        out << "******************************* HEADER *******************************" << endl;
        header_->Print(out);
    }
    if(title_ != NULL)
    {
        out << "******************************** TITLE *******************************" << endl;
        title_->Print(out);
    }
    if(compound_ != NULL)
    {
        out << "****************************** COMPOUND ******************************" << endl;
        compound_->Print(out);
    }
    if(number_of_models_ != NULL)
    {
        out << "************************** NUMBER OF MODELS **************************" << endl;
        number_of_models_->Print(out);
    }
    if(model_type_ != NULL)
    {
        out << "***************************** MODEL TYPE *****************************" << endl;
        model_type_->Print(out);
    }
    if(residues_sequence_ != NULL)
    {
        out << "************************** RESIDUE SEQUENCE **************************" << endl;
        residues_sequence_->Print(out);
    }
    if(residue_modification_ != NULL)
    {
        out << "************************ RESIDUE MODIFICATION ************************" << endl;
        residue_modification_->Print(out);
    }
    if(heterogens_ != NULL)
    {
        out << "***************************** HETEROGEN ******************************" << endl;
        heterogens_->Print(out);
    }
    if(heterogens_name_ != NULL)
    {
        out << "*************************** HETEROGEN NAME ***************************" << endl;
        heterogens_name_->Print(out);
    }
    if(heterogen_synonyms_ != NULL)
    {
        out << "************************** HETEROGEN SYNONYM *************************" << endl;
        heterogen_synonyms_->Print(out);
    }
    if(formulas_ != NULL)
    {
        out << "******************************* FORMULA ******************************" << endl;
        formulas_->Print(out);
    }
    if(helixes_ != NULL)
    {
        out << "******************************** HELIX *******************************" << endl;
        helixes_->Print(out);
    }
    if(sheets_ != NULL)
    {
        out << "******************************** SHEET *******************************" << endl;
        sheets_->Print(out);
    }
    if(disulfide_bonds_ != NULL)
    {
        out << "*************************** DISULFIDE BOND ***************************" << endl;
        disulfide_bonds_->Print(out);
    }
    if(links_ != NULL)
    {
        out << "******************************** LINK ********************************" << endl;
        links_->Print(out);
    }
    if(sites_ != NULL)
    {
        out << "******************************** SITE ********************************" << endl;
        sites_->Print(out);
    }
    if(crystallography_ != NULL)
    {
        out << "************************** CRYSTALLOGRAPHIC **************************" << endl;
        crystallography_->Print(out);
    }
    if(origins_ != NULL)
    {
        out << "******************************* ORIGIN *******************************" << endl;
        origins_->Print(out);
    }
    if(scales_ != NULL)
    {
        out << "******************************** SCALE *******************************" << endl;
        scales_->Print(out);
    }
    if(matrices_ != NULL)
    {
        out << "******************************* MATRIX *******************************" << endl;
        matrices_->Print(out);
    }
    if(models_ != NULL)
    {
        out << "******************************* MODEL ********************************" << endl;
        models_->Print(out);
    }
    if(connectivities_ != NULL)
    {
        out << "******************************* CONNECT ******************************" << endl;
        connectivities_->Print(out);
    }
}
