#include <cctype>
#include <ctime>
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessor.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "includes/CodeUtils/logging.hpp"


using PdbPreprocessorSpace::PdbPreprocessor;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessor::PdbPreprocessor(PdbFileSpace::PdbFile &pdbFile) : pdbFile_ (pdbFile)
{
    this->Preprocess();
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbPreprocessor::PdbPreprocessorDisulfideBondVector PdbPreprocessor::GetDisulfideBonds(){
    return disulfide_bonds_;
}
PdbPreprocessor::PdbPreprocessorChainTerminationVector PdbPreprocessor::GetChainTerminations(){
    return chain_terminations_;
}
PdbPreprocessor::PdbPreprocessorHistidineMappingVector PdbPreprocessor::GetHistidineMappings(){
    return histidine_mappings_;
}
PdbPreprocessor::PdbPreprocessorMissingResidueVector PdbPreprocessor::GetMissingResidues(){
    return missing_residues_;
}
PdbPreprocessor::PdbPreprocessorUnrecognizedResidueVector PdbPreprocessor::GetUnrecognizedResidues(){
    return unrecognized_residues_;
}
PdbPreprocessor::PdbPreprocessorUnrecognizedHeavyAtomVector PdbPreprocessor::GetUnrecognizedHeavyAtoms(){
    return unrecognized_heavy_atoms_;
}
PdbPreprocessor::PdbPreprocessorReplacedHydrogenVector PdbPreprocessor::GetReplacedHydrogens(){
    return replaced_hydrogens_;
}
PdbPreprocessor::PdbPreprocessorAlternateResidueMap PdbPreprocessor::GetAlternateResidueMap(){
    return alternate_residue_map_;
}
PdbPreprocessor::PdbPreprocessorToBeDeletedAtomVector PdbPreprocessor::GetToBeDeletedAtoms(){
    return to_be_deleted_atoms_;
}
PdbPreprocessor::PdbPreprocessorToBeDeletedResidueVector PdbPreprocessor::GetToBeDeletedResidues(){
    return to_be_deleted_residues_;
}
PdbPreprocessor::PdbPreprocessorResidueInfoMap PdbPreprocessor::GetResidueInfoMap(){
    return residue_info_map_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

void PdbPreprocessor::AddDisulfideBond(PdbPreprocessorDisulfideBond *disulfide_bond)
{
    disulfide_bonds_.push_back(disulfide_bond);
}
void PdbPreprocessor::AddChainTermination(PdbPreprocessorChainTermination *chain_termination)
{
    chain_terminations_.push_back(chain_termination);
}
void PdbPreprocessor::AddHistidineMapping(PdbPreprocessorHistidineMapping *histidine_mapping)
{
    histidine_mappings_.push_back(histidine_mapping);
}
void PdbPreprocessor::AddMissingResidue(PdbPreprocessorMissingResidue *missing_residue)
{
    missing_residues_.push_back(missing_residue);
}
void PdbPreprocessor::AddUnrecognizedResidue(PdbPreprocessorUnrecognizedResidue *unrecognized_residue)
{
    unrecognized_residues_.push_back(unrecognized_residue);
}
void PdbPreprocessor::AddUnrecognizedHeavyAtom(PdbPreprocessorUnrecognizedHeavyAtom *unrecognized_heavy_atom)
{
    unrecognized_heavy_atoms_.push_back(unrecognized_heavy_atom);
}
void PdbPreprocessor::AddReplacedHydrogen(PdbPreprocessorReplacedHydrogen *replaced_hydrogen)
{
    replaced_hydrogens_.push_back(replaced_hydrogen);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

gmml::ResidueNameMap PdbPreprocessor::GetUnrecognizedResidueNamesMap(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, gmml::ResidueNameMap dataset_residue_names)
{
    gmml::ResidueNameMap unrecognized_residue_names;
    std::pair<std::string, std::string> residue_name_flag;
    std::string pdb_residue_name;
    std::string pdb_residue_position_flag;
    for(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        residue_name_flag = *it;
        pdb_residue_name = residue_name_flag.first;
        pdb_residue_position_flag = residue_name_flag.second;
        if(pdb_residue_name.compare("HIS") != 0 )
        {
            if(dataset_residue_names.find(pdb_residue_name) != dataset_residue_names.end())
            {
            }
            else
            {
                if(pdb_residue_position_flag.compare("S") == 0)
                {
                    std::stringstream ss;
                    ss << "N" << pdb_residue_name;
                    std::string n_terminal_pdb_residue_name = ss.str();
                    if(dataset_residue_names.find(n_terminal_pdb_residue_name) == dataset_residue_names.end())
                        unrecognized_residue_names[pdb_residue_name] = pdb_residue_name;
                }
                else if(pdb_residue_position_flag.compare("E") == 0)
                {
                    std::stringstream ss;
                    ss << "C" << pdb_residue_name;
                    std::string c_terminal_pdb_residue_name = ss.str();
                    if(dataset_residue_names.find(c_terminal_pdb_residue_name) == dataset_residue_names.end())
                        unrecognized_residue_names[pdb_residue_name] = pdb_residue_name;
                }
                else
                {
                    unrecognized_residue_names[pdb_residue_name] = pdb_residue_name;
                }
            }
        }
    }
//    std::cout << "HIS residue(s) found" << std::endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "HIS residue(s) found" ); // wtf?
    return unrecognized_residue_names;
}

gmml::ResidueNameMap PdbPreprocessor::GetRecognizedResidueNamesMap(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, gmml::ResidueNameMap dataset_residue_names)
{
    gmml::ResidueNameMap recognized_residue_names = gmml::ResidueNameMap();
    std::pair<std::string, std::string> residue_name_flag;
    std::string pdb_residue_name;
    std::string pdb_residue_position_flag;
    for(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        residue_name_flag = *it;
        pdb_residue_name = residue_name_flag.first;
        pdb_residue_position_flag = residue_name_flag.second;
        if(pdb_residue_name.compare("HIS") != 0)
        {
            if(dataset_residue_names.find(pdb_residue_name) != dataset_residue_names.end())
            {
                recognized_residue_names[pdb_residue_name] = pdb_residue_name;
            }
            else
            {
                if(pdb_residue_position_flag.compare("S") == 0)
                {
                    std::stringstream ss;
                    ss << "N" << pdb_residue_name;
                    std::string n_terminal_pdb_residue_name = ss.str();
                    if(dataset_residue_names.find(n_terminal_pdb_residue_name) != dataset_residue_names.end())
                        recognized_residue_names[pdb_residue_name] = pdb_residue_name;
                }
                else if(pdb_residue_position_flag.compare("E") == 0)
                {
                    std::stringstream ss;
                    ss << "C" << pdb_residue_name;
                    std::string c_terminal_pdb_residue_name = ss.str();
                    if(dataset_residue_names.find(c_terminal_pdb_residue_name) != dataset_residue_names.end())
                        recognized_residue_names[pdb_residue_name] = pdb_residue_name;
                }
            }
        }
        else
        {
            recognized_residue_names[pdb_residue_name] = pdb_residue_name;
        }
    }
//    std::cout << "HIS residue(s) found" << std::endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "HIS residue(s) found" );
    return recognized_residue_names;
}

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetUnrecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, gmml::ResidueNameMap unrecognized_residue_names)
{
    PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* pdb_residue = *it;
        std::string pdb_residue_name = pdb_residue->GetResidueName();
        if(unrecognized_residue_names.find(pdb_residue_name) != unrecognized_residue_names.end())
            unrecognized_residues.push_back(pdb_residue);
    }
    return unrecognized_residues;
}

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetRecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, gmml::ResidueNameMap recognized_residue_names)
{
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* pdb_residue = *it;
        std::string pdb_residue_name = pdb_residue->GetResidueName();
        if(recognized_residue_names.find(pdb_residue_name) != recognized_residue_names.end())
            recognized_residues.push_back(pdb_residue);
    }
    return recognized_residues;
}

bool PdbPreprocessor::ExtractUnrecognizedResidues()
{
    // Advanced version
    std::map<std::string, std::string> dataset_residue_names = this->GetParameters().GetAllResidueNameMap();
    std::vector<std::pair<std::string, std::string> > pdb_residue_names = this->GetPdbFile().GetAllResidueNames();
    std::map<std::string, std::string> unrecognized_residue_names = GetUnrecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);

    std::vector<PdbFileSpace::PdbResidue*> pdb_residues = this->GetPdbFile().GetAllResidues();
    std::vector<PdbFileSpace::PdbResidue*> unrecognized_residues = GetUnrecognizedResidues(pdb_residues, unrecognized_residue_names);

    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues_atom_card = this->GetPdbFile().GetAllResiduesFromAtomSection();
    PdbPreprocessorChainIdResidueMap all_chain_map_residue;
    PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
    PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues_atom_card.begin(); it != pdb_residues_atom_card.end(); it++)
    {
        PdbFileSpace::PdbResidue* residue = *it;
        char chain_id = residue->GetResidueChainId();

        std::stringstream ss;
        ss << chain_id;
        all_chain_map_residue[ss.str()].push_back(residue);
    }
    gmml::ResidueNameMap amino_lib_residue_names_map = this->GetParameters().GetAllResidueNameMap();
    PdbPreprocessorChainIdResidueMap chain_map_residue;
    PdbFileSpace::PdbFile::PdbResidueVector residues;
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues_in_chain;
    for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
    {
        int internal_amino_acid_chain_counter = 0;
        chain_map_residue = PdbPreprocessorChainIdResidueMap();
        std::string chain_id = (*it).first;
        residues = (*it).second;

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            PdbFileSpace::PdbResidue* residue = (*it1);
            if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
            {
                std::stringstream ss;
                ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                chain_map_residue[ss.str()].push_back(residue);
                std::stringstream sss(chain_id);
                char c_id;
                sss >> c_id;
                chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
            }
            else
            {
                std::stringstream ss;
                ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                chain_map_residue[ss.str()].push_back(residue);
                internal_amino_acid_chain_counter++;
            }
        }

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = unrecognized_residues.begin(); it1 != unrecognized_residues.end(); it1++)
        {
            PdbFileSpace::PdbResidue* pdb_residue = (*it1);
            PdbPreprocessorUnrecognizedResidue* unrecognized_residue =
                    new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(), pdb_residue->GetResidueSequenceNumber(),
                                                           pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation(), false);
            if(chain_map_residue.size() <= 2)
                unrecognized_residue->SetMiddleOfChain(false);
            else
            {
                for(PdbPreprocessorChainIdResidueMap::iterator it2 = chain_map_residue.begin(); it2 != chain_map_residue.end(); it2++)
                {
                    std::string key = (*it2).first;
                    pdb_residues_in_chain = (*it2).second;
                    unsigned int dist = distance(chain_map_residue.begin(), it2);
                    if(key.substr(0, 2).compare("NA") == 0 && dist != 0 && dist != chain_map_residue.size() - 1)
                    {
                        if(find(pdb_residues_in_chain.begin(), pdb_residues_in_chain.end(), pdb_residue) != pdb_residues_in_chain.end())
                        {
                            unrecognized_residue->SetMiddleOfChain(true);
                        }
                    }
                    else
                        unrecognized_residue->SetMiddleOfChain(false);
                }
            }
            this->AddUnrecognizedResidue(unrecognized_residue);
        }
    }
    return true;
}

void PdbPreprocessor::RemoveUnrecognizedResiduesWithTheGivenModelNumber(PdbPreprocessorUnrecognizedResidueVector unrecognized_residues, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorUnrecognizedResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
    {
        PdbPreprocessorUnrecognizedResidue* unrecognized_residue = (*it);
        PdbFileSpace::PdbResidue* pdb_residue =
                new PdbFileSpace::PdbResidue(unrecognized_residue->GetResidueName(), unrecognized_residue->GetResidueChainId(),
                               unrecognized_residue->GetResidueSequenceNumber(), unrecognized_residue->GetResidueInsertionCode(), unrecognized_residue->GetResidueAlternateLocation());
//        this->GetPdbFile().DeleteResidueWithTheGivenModelNumber(pdb_residue, model_number);
        to_be_deleted_residues_.push_back(pdb_residue);
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(model_number);
}

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetAllCYSResidues(PdbFileSpace::PdbFile::PdbResidueVector pdbResidues)
{
    PdbFileSpace::PdbFile::PdbResidueVector allCysResidues;
    for(auto &pdbResidue : pdbResidues)
    {
        if(pdbResidue->GetResidueName() == "CYS")
        {
            allCysResidues.push_back(pdbResidue);
        }
    }
    return allCysResidues;
}

double PdbPreprocessor::GetDistanceofCYS(PdbFileSpace::PdbResidue *first_residue, PdbFileSpace::PdbResidue *second_residue,
                                         PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atom_map, int &first_sulfur_atom_serial_number,
                                         int &second_sulfur_atom_serial_number)
{
    double distance = 0.0;
    PdbFileSpace::PdbAtomCard* first_residue_sulfur_atom = this->GetPdbFile().GetAtomOfResidueByName(first_residue, "SG", residue_atom_map);
    PdbFileSpace::PdbAtomCard* second_residue_sulfur_atom = this->GetPdbFile().GetAtomOfResidueByName(second_residue, "SG", residue_atom_map);
    first_sulfur_atom_serial_number = first_residue_sulfur_atom->GetAtomSerialNumber();
    second_sulfur_atom_serial_number = second_residue_sulfur_atom->GetAtomSerialNumber();
    if(first_residue_sulfur_atom != NULL && second_residue_sulfur_atom != NULL)
        distance = first_residue_sulfur_atom->GetAtomOrthogonalCoordinate().Distance(second_residue_sulfur_atom->GetAtomOrthogonalCoordinate());
    return distance;
}

bool PdbPreprocessor::ExtractCYSResidues()
{
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResiduesFromAtomSection();
    PdbFileSpace::PdbFile::PdbResidueVector cys_residues = GetAllCYSResidues(pdb_residues);
    PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atom_map = this->GetPdbFile().GetAllAtomsOfResidues();

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = cys_residues.begin(); it != cys_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* first_residue = (*it);
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = it+1 ; it1 != cys_residues.end(); it1++)
        {
            PdbFileSpace::PdbResidue* second_residue = (*it1);
            int first_sulfur_atom_serial_number;
            int second_sulfur_atom_serial_number;
            double distance = GetDistanceofCYS(first_residue, second_residue, residue_atom_map, first_sulfur_atom_serial_number,
                                               second_sulfur_atom_serial_number);
            if (distance < gmml::dSulfurCutoff)
            {
                PdbPreprocessorDisulfideBond* disulfide_bond =
                        new PdbPreprocessorDisulfideBond(first_residue->GetResidueChainId(), second_residue->GetResidueChainId(),
                                                         first_residue->GetResidueSequenceNumber(), second_residue->GetResidueSequenceNumber(),
                                                         distance, true, first_residue->GetResidueInsertionCode(), second_residue->GetResidueInsertionCode(),
                                                         first_residue->GetResidueAlternateLocation(), second_residue->GetResidueAlternateLocation(),
                                                         first_sulfur_atom_serial_number, second_sulfur_atom_serial_number);
                this->AddDisulfideBond(disulfide_bond);
            }
        }
    }
    return true;
}

void PdbPreprocessor::UpdateCYSResiduesWithTheGivenModelNumber(PdbPreprocessorDisulfideBondVector disulfide_bonds, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues;
    for(PdbPreprocessorDisulfideBondVector::iterator it = disulfide_bonds.begin(); it != disulfide_bonds.end(); it++)
    {
        PdbPreprocessorDisulfideBond* disulfide_bond = (*it);
        if(disulfide_bond->GetIsBonded())
        {
            PdbFileSpace::PdbAtomCard* pdb_atom_1 =
                    new PdbFileSpace::PdbAtomCard(disulfide_bond->GetResidueChainId1(), "HG", "CYS",
                                disulfide_bond->GetResidueSequenceNumber1(), disulfide_bond->GetResidueInsertionCode1(), disulfide_bond->GetResidueAlternateLocation1());
//            this->GetPdbFile().DeleteAtomWithTheGivenModelNumber(pdb_atom_1, model_number);
            to_be_deleted_atoms_.push_back(pdb_atom_1);
            PdbFileSpace::PdbAtomCard* pdb_atom_2 =
                    new PdbFileSpace::PdbAtomCard(disulfide_bond->GetResidueChainId2(), "HG", "CYS",
                                disulfide_bond->GetResidueSequenceNumber2(), disulfide_bond->GetResidueInsertionCode2(), disulfide_bond->GetResidueAlternateLocation2());
//            this->GetPdbFile().DeleteAtomWithTheGivenModelNumber(pdb_atom_2, model_number);
            to_be_deleted_atoms_.push_back(pdb_atom_2);
            std::string target_key1;
            std::stringstream ss_1;
            ss_1 << "CYS" << "_" << disulfide_bond->GetResidueChainId1() << "_" << disulfide_bond->GetResidueSequenceNumber1() << "_" << disulfide_bond->GetResidueInsertionCode1()
                 << "_" << disulfide_bond->GetResidueAlternateLocation1();
            target_key1 = ss_1.str();
            std::string target_key2;
            std::stringstream ss_2;
            ss_2 << "CYS" << "_" << disulfide_bond->GetResidueChainId2() << "_" << disulfide_bond->GetResidueSequenceNumber2() << "_" << disulfide_bond->GetResidueInsertionCode2()
                 << "_" << disulfide_bond->GetResidueAlternateLocation2();
            target_key2 = ss_2.str();
            pdb_residues = this->GetPdbFile().GetAllResidues();
            for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = pdb_residues.begin(); it1 != pdb_residues.end(); it1++)
            {
                PdbFileSpace::PdbResidue* pdb_residue = (*it1);

                std::string residue_name = pdb_residue->GetResidueName();
                if((residue_name).compare("CYS") == 0)
                {
                    std::stringstream ss;
                    std::string pdb_residue_key;
                    ss << residue_name << "_" << pdb_residue->GetResidueChainId() << "_" << pdb_residue->GetResidueSequenceNumber() << "_" << pdb_residue->GetResidueInsertionCode()
                       << "_" << pdb_residue->GetResidueAlternateLocation();
                    pdb_residue_key = ss.str();
                    if(pdb_residue_key.compare(target_key1) == 0 || pdb_residue_key.compare(target_key2) == 0)
                    {
                        this->GetPdbFile().UpdateResidueNameWithTheGivenModelNumber(pdb_residue, "CYX", model_number);
                    }
                }
            }
            PdbFileSpace::PdbConnectSection* pdb_connect_card;
            if(this->GetPdbFile().GetConnectivities() == NULL)
            {
                pdb_connect_card = new PdbFileSpace::PdbConnectSection();
            }
            else
            {
                pdb_connect_card = this->GetPdbFile().GetConnectivities();
            }
            PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap bonded_atoms_serial_numbers_map = pdb_connect_card->GetBondedAtomsSerialNumbers();
            if(bonded_atoms_serial_numbers_map.find(disulfide_bond->GetSulfurAtomSerialNumber1()) != bonded_atoms_serial_numbers_map.end())
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()].push_back(disulfide_bond->GetSulfurAtomSerialNumber2());
            else
            {
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()] = std::vector<int>();
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()].push_back(disulfide_bond->GetSulfurAtomSerialNumber2());
            }
            if(bonded_atoms_serial_numbers_map.find(disulfide_bond->GetSulfurAtomSerialNumber2()) != bonded_atoms_serial_numbers_map.end())
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()].push_back(disulfide_bond->GetSulfurAtomSerialNumber1());
            else
            {
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()] = std::vector<int>();
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()].push_back(disulfide_bond->GetSulfurAtomSerialNumber1());
            }

            pdb_connect_card->SetBondedAtomsSerialNumbers(bonded_atoms_serial_numbers_map);
            this->GetPdbFile().SetConnectivities(pdb_connect_card);
        }
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(model_number);
}

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetAllHISResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues)
{
    PdbFileSpace::PdbFile::PdbResidueVector all_his_residues;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* pdb_residue = (*it);
        std::string pdb_residue_name = pdb_residue->GetResidueName();
        if((pdb_residue_name).compare("HIS") == 0)
        {
            all_his_residues.push_back(pdb_residue);
        }
    }
    return all_his_residues;
}

bool PdbPreprocessor::ExtractHISResidues()
{
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResiduesFromAtomSection();
    PdbFileSpace::PdbFile::PdbResidueVector his_residues = GetAllHISResidues(pdb_residues);
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = his_residues.begin(); it != his_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* his_residue = (*it);

        // gmml::HIE residue
        if(this->GetPdbFile().GetAtomOfResidueByName(his_residue, "HE2") != NULL && this->GetPdbFile().GetAtomOfResidueByName(his_residue, "HD1") == NULL)
        {
            PdbPreprocessorHistidineMapping* histidine_mapping =
                    new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), gmml::HIE,
                                                        his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
            this->AddHistidineMapping(histidine_mapping);
        }
        // gmml::HID residue
        else if(this->GetPdbFile().GetAtomOfResidueByName(his_residue, "HE2") == NULL && this->GetPdbFile().GetAtomOfResidueByName(his_residue, "HD1") != NULL)
        {
            PdbPreprocessorHistidineMapping* histidine_mapping =
                    new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), gmml::HID,
                                                        his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
            this->AddHistidineMapping(histidine_mapping);
        }
        // gmml::HIP residue
        else if(this->GetPdbFile().GetAtomOfResidueByName(his_residue, "HE2") != NULL && this->GetPdbFile().GetAtomOfResidueByName(his_residue, "HD1") != NULL)
        {
            PdbPreprocessorHistidineMapping* histidine_mapping =
                    new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), gmml::HIP,
                                                        his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
            this->AddHistidineMapping(histidine_mapping);
        }
        else if(this->GetPdbFile().GetAtomOfResidueByName(his_residue, "HE2") == NULL && this->GetPdbFile().GetAtomOfResidueByName(his_residue, "HD1") == NULL)
        {
            PdbPreprocessorHistidineMapping* histidine_mapping =
                    new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), gmml::HIE,
                                                        his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
            this->AddHistidineMapping(histidine_mapping);
        }
    }
    return true;
}

void PdbPreprocessor::UpdateHISMappingWithTheGivenNumber(PdbPreprocessorHistidineMappingVector histidine_mappings, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResidues();
    for(PdbPreprocessorHistidineMappingVector::iterator it = histidine_mappings.begin(); it != histidine_mappings.end(); it++)
    {
        PdbPreprocessorHistidineMapping* histidine_mapping = (*it);
        std::string target_key;
        std::stringstream ss;
        ss << "HIS" << "_" << histidine_mapping->GetResidueChainId() << "_" << histidine_mapping->GetResidueSequenceNumber()
           << "_" << histidine_mapping->GetResidueInsertionCode() << "_" << histidine_mapping->GetResidueAlternateLocation();
        target_key = ss.str();
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = pdb_residues.begin(); it1 != pdb_residues.end(); it1++)
        {
            PdbFileSpace::PdbResidue* pdb_residue = (*it1);

            std::string residue_name = pdb_residue->GetResidueName();
            if((residue_name).compare("HIS") == 0)
            {
                std::stringstream ss1;
                std::string pdb_residue_key;
                ss1 << residue_name << "_" << pdb_residue->GetResidueChainId() << "_" << pdb_residue->GetResidueSequenceNumber()
                    << "_" << pdb_residue->GetResidueInsertionCode() << "_" << pdb_residue->GetResidueAlternateLocation();
                pdb_residue_key = ss1.str();
                if(pdb_residue_key.compare(target_key) == 0)
                {
                    // gmml::HIE residue
                    PdbFileSpace::PdbAtomCard* HE2 = this->GetPdbFile().GetAtomOfResidueByName(pdb_residue, "HE2");
                    PdbFileSpace::PdbAtomCard* HD1 = this->GetPdbFile().GetAtomOfResidueByName(pdb_residue, "HD1");
                    if(HE2 != NULL && HD1 == NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == gmml::HID)
                        {
                            // Delete HE2
//                            this->GetPdbFile().DeleteAtomWithTheGivenModelNumber(HE2, model_number);
                            to_be_deleted_atoms_.push_back(HE2);
                        }
                    }
                    // gmml::HID residue
                    if(HE2 == NULL && HD1 != NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == gmml::HIE)
                        {
                            // Delete HD1
//                            this->GetPdbFile().DeleteAtomWithTheGivenModelNumber(HD1, model_number);
                            to_be_deleted_atoms_.push_back(HD1);
                        }

                    }
                    // gmml::HIP residue
                    if(HE2 != NULL && HD1 != NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == gmml::HIE)
                        {
                            // Delete HD1
//                            this->GetPdbFile().DeleteAtomWithTheGivenModelNumber(HD1, model_number);
                            to_be_deleted_atoms_.push_back(HD1);
                        }
                        if(histidine_mapping->GetSelectedMapping() == gmml::HID)
                        {
                            // Delete HE2
//                            this->GetPdbFile().DeleteAtomWithTheGivenModelNumber(HE2, model_number);
                            to_be_deleted_atoms_.push_back(HE2);
                        }
                    }
                    this->GetPdbFile().UpdateResidueNameWithTheGivenModelNumber(pdb_residue, histidine_mapping->GetStringFormatOfSelectedMapping(), model_number);
                }
            }
        }
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(model_number);
}

PdbFileSpace::PdbFile::PdbAtomCardVector PdbPreprocessor::GetUnknownHeavyAtomsOfResidue(PdbFileSpace::PdbFile::PdbAtomCardVector pdb_atoms, std::vector<std::string> dataset_atom_names_of_residue)
{
    PdbFileSpace::PdbFile::PdbAtomCardVector unknown_heavy_atoms_of_residue;
    for(PdbFileSpace::PdbFile::PdbAtomCardVector::iterator it = pdb_atoms.begin(); it != pdb_atoms.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* pdb_atom = *it;
        std::string pdb_atom_name = pdb_atom->GetAtomName();
        if(!(pdb_atom_name.substr(0,1).compare("H") == 0 ||
             (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(gmml::ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(!(find(dataset_atom_names_of_residue.begin(), dataset_atom_names_of_residue.end(), pdb_atom_name) != dataset_atom_names_of_residue.end()))
            {
                unknown_heavy_atoms_of_residue.push_back(pdb_atom);
            }
        }
    }
    return unknown_heavy_atoms_of_residue;
}

bool PdbPreprocessor::ExtractUnknownHeavyAtoms()
{
    gmml::ResidueNameMap dataset_residue_names = this->GetParameters().GetAllResidueNameMap();
    std::vector<std::pair<std::string, std::string> > pdb_residue_names = this->GetPdbFile().GetAllResidueNames(); // It will be the residue name and then either S for start, E for end residue of each model in the PDB file, "" for nothing. Not great.
    gmml::ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);


    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

    PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atom_map = this->GetPdbFile().GetAllAtomsOfResidues();
    std::map<std::string, std::vector<std::string> > dataset_residue_atom_map = this->GetParameters().GetResidueNamesToTheirAtomNamesMap();
    std::vector<std::string> dataset_atom_names_of_residue = std::vector<std::string>();
    std::vector<std::string> dataset_atom_names_of_head_residue = std::vector<std::string>();
    std::vector<std::string> dataset_atom_names_of_tail_residue = std::vector<std::string>();
    PdbFileSpace::PdbFile::PdbAtomCardVector unknown_heavy_atoms = PdbFileSpace::PdbFile::PdbAtomCardVector();
    PdbFileSpace::PdbFile::PdbAtomCardVector atoms_of_residue = PdbFileSpace::PdbFile::PdbAtomCardVector();
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* recognized_residue = *it;
        std::string residue_name = recognized_residue->GetResidueName();
        if(residue_name.compare("HIS") != 0)
        {
            char chain_id = recognized_residue->GetResidueChainId();
            int sequence_number = recognized_residue->GetResidueSequenceNumber();
            char insertion_code = recognized_residue->GetResidueInsertionCode();
            char alternate_location = recognized_residue->GetResidueAlternateLocation();

            std::stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            std::string key = ss.str();
            atoms_of_residue = *(residue_atom_map[key]);
            dataset_atom_names_of_residue = dataset_residue_atom_map[residue_name];

            std::pair<std::string, std::string> residue_sflag_pair = std::make_pair(residue_name, "S");
            std::pair<std::string, std::string> residue_eflag_pair = std::make_pair(residue_name, "E");
            if(find(pdb_residue_names.begin(), pdb_residue_names.end(), residue_sflag_pair) != pdb_residue_names.end())
            {
                std::stringstream ss1;
                ss1 << "N" << residue_name;
                dataset_atom_names_of_head_residue = dataset_residue_atom_map[ss1.str()];
                for(std::vector<std::string>::iterator it1 = dataset_atom_names_of_head_residue.begin(); it1 != dataset_atom_names_of_head_residue.end(); it1++)
                {
                    std::string head_residue = (*it1);
                    dataset_atom_names_of_residue.push_back(head_residue);
                }
            }
            else if(find(pdb_residue_names.begin(), pdb_residue_names.end(), residue_eflag_pair) != pdb_residue_names.end())
            {
                std::stringstream ss1;
                ss1 << "C" << residue_name;
                dataset_atom_names_of_tail_residue = dataset_residue_atom_map[ss1.str()];
                for(std::vector<std::string>::iterator it2 = dataset_atom_names_of_tail_residue.begin(); it2 != dataset_atom_names_of_tail_residue.end(); it2++)
                {
                    std::string head_residue = (*it2);
                    dataset_atom_names_of_residue.push_back(head_residue);
                }
            }

            unknown_heavy_atoms = GetUnknownHeavyAtomsOfResidue(atoms_of_residue, dataset_atom_names_of_residue);

            for(PdbFileSpace::PdbFile::PdbAtomCardVector::iterator it1 = unknown_heavy_atoms.begin(); it1 != unknown_heavy_atoms.end(); it1++)
            {
                PdbFileSpace::PdbAtomCard* heavy_atom = (*it1);
                PdbPreprocessorUnrecognizedHeavyAtom* unknown_heavy_atom =
                        new PdbPreprocessorUnrecognizedHeavyAtom(heavy_atom->GetAtomChainId(), heavy_atom->GetAtomSerialNumber(),
                                                                 heavy_atom->GetAtomName(), heavy_atom->GetAtomResidueName(),
                                                                 heavy_atom->GetAtomResidueSequenceNumber(), heavy_atom->GetAtomInsertionCode(), heavy_atom->GetAtomAlternateLocation());
                this->AddUnrecognizedHeavyAtom(unknown_heavy_atom);
            }
        }
    }
    return true;
}

void PdbPreprocessor::RemoveResiduesOfUnknownHeavyAtomsWithTheGivenModelNumber(PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    std::vector<std::string> removed_keys = std::vector<std::string>();
    for(PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unknown_heavy_atoms.begin(); it != unknown_heavy_atoms.end(); it++)
    {
        PdbPreprocessorUnrecognizedHeavyAtom* unknown_heavy_atom = (*it);
        std::stringstream ss;
        ss << unknown_heavy_atom->GetResidueInsertionCode() << "_" << unknown_heavy_atom->GetResidueChainId() << "_" << unknown_heavy_atom->GetResidueSequenceNumber()
           << "_" << unknown_heavy_atom->GetResidueInsertionCode() << "_" << unknown_heavy_atom->GetResidueAlternateLocation();
        std::string residue_key = ss.str();
        if(distance(removed_keys.begin(), find(removed_keys.begin(), removed_keys.end(), residue_key)) < 0 ||
                distance(removed_keys.begin(), find(removed_keys.begin(), removed_keys.end(), residue_key)) >= (int)removed_keys.size())
        {
            PdbFileSpace::PdbResidue* pdb_residue = new PdbFileSpace::PdbResidue(unknown_heavy_atom->GetResidueName(),unknown_heavy_atom->GetResidueChainId(), unknown_heavy_atom->GetResidueSequenceNumber(),
                                                     unknown_heavy_atom->GetResidueInsertionCode(), unknown_heavy_atom->GetResidueAlternateLocation());
//            this->GetPdbFile().DeleteResidueWithTheGivenModelNumber(pdb_residue , model_number);
            to_be_deleted_residues_.push_back(pdb_residue);
            removed_keys.push_back(residue_key);
        }
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(model_number);
}

PdbFileSpace::PdbFile::PdbAtomCardVector PdbPreprocessor::GetRemovedHydrogensOfResidue(PdbFileSpace::PdbFile::PdbAtomCardVector pdb_atoms, std::vector<std::string> dataset_atom_names_of_residue)
{
    PdbFileSpace::PdbFile::PdbAtomCardVector removed_hydrogens_of_residue;
    for(PdbFileSpace::PdbFile::PdbAtomCardVector::iterator it = pdb_atoms.begin(); it != pdb_atoms.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* pdb_atom = *it;
        std::string pdb_atom_name = pdb_atom->GetAtomName();
        if((pdb_atom_name.substr(0,1).compare("H") == 0 ||
            (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(gmml::ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(find(dataset_atom_names_of_residue.begin(), dataset_atom_names_of_residue.end(), pdb_atom_name) == dataset_atom_names_of_residue.end())
            {
                removed_hydrogens_of_residue.push_back(pdb_atom);
                std::stringstream ss;
                ss << pdb_atom->GetAtomName() << "_" << pdb_atom->GetAtomResidueSequenceNumber() << "_" << pdb_atom->GetAtomResidueName() << "_" << pdb_atom->GetAtomChainId();
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "Removing hydrogen: " + ss.str());
            }
        }
    }
    return removed_hydrogens_of_residue;
}

bool PdbPreprocessor::ExtractRemovedHydrogens()
{
    gmml::ResidueNameMap dataset_residue_names = this->GetParameters().GetAllResidueNameMap();
    PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = this->GetPdbFile().GetAllResidueNames();
    gmml::ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

    PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atom_map = this->GetPdbFile().GetAllAtomsOfResidues();
    PdbFileSpace::PdbFile::PdbAtomCardVector atoms_of_residue;
    gmml::ResidueNameAtomNamesMap dataset_residue_atom_map = this->GetParameters().GetResidueNamesToTheirAtomNamesMap();
    PdbFileSpace::PdbFile::PdbAtomCardVector removed_hydrogens;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* recognized_residue = *it;
        if(recognized_residue->GetResidueName().compare("HIS") != 0)
        {
            std::string residue_name = recognized_residue->GetResidueName();
            char chain_id = recognized_residue->GetResidueChainId();
            int sequence_number = recognized_residue->GetResidueSequenceNumber();
            char insertion_code = recognized_residue->GetResidueInsertionCode();
            char alternate_location = recognized_residue->GetResidueAlternateLocation();
            std::stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            std::string key = ss.str();
            atoms_of_residue = *(residue_atom_map[key]);
            std::vector<std::string> dataset_atom_names_of_residue = dataset_residue_atom_map[residue_name];
            // OG JUST DUCT TAPING AT THIS POINT. Sometimes there will be atoms coming in like H1, H2, H3 on N terminus, we want to keep those.
            std::vector<std::string> dataset_atom_names_of_Nresidue = dataset_residue_atom_map[("N" + residue_name)];
            std::vector<std::string> dataset_atom_names_of_Cresidue = dataset_residue_atom_map[("C" + residue_name)];
            dataset_atom_names_of_residue.insert(dataset_atom_names_of_residue.end(), dataset_atom_names_of_Nresidue.begin(), dataset_atom_names_of_Nresidue.end());
            dataset_atom_names_of_residue.insert(dataset_atom_names_of_residue.end(), dataset_atom_names_of_Cresidue.begin(), dataset_atom_names_of_Cresidue.end());
            // OG END DUCT TAPE.
            removed_hydrogens = GetRemovedHydrogensOfResidue(atoms_of_residue, dataset_atom_names_of_residue);
            for(PdbFileSpace::PdbFile::PdbAtomCardVector::iterator it1 = removed_hydrogens.begin(); it1 != removed_hydrogens.end(); it1++)
            {
                PdbFileSpace::PdbAtomCard* removed_hydrogen = (*it1);

                gmml::log(__LINE__, __FILE__, gmml::INF, "Removing hydrogen:" + removed_hydrogen->GetAtomName());
                PdbPreprocessorReplacedHydrogen* removed_hydrogen_atom =
                        new PdbPreprocessorReplacedHydrogen(removed_hydrogen->GetAtomChainId(), removed_hydrogen->GetAtomSerialNumber(), removed_hydrogen->GetAtomName(),
                                                            removed_hydrogen->GetAtomResidueName(), removed_hydrogen->GetAtomResidueSequenceNumber(),
                                                            removed_hydrogen->GetAtomInsertionCode(), removed_hydrogen->GetAtomAlternateLocation());
                this->AddReplacedHydrogen(removed_hydrogen_atom);
            }
        }
    }
    return true;
}

void PdbPreprocessor::RemoveRemovedHydrogensWithTheGivenModelNumber(PdbPreprocessorReplacedHydrogenVector replaced_hydrogens, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorReplacedHydrogenVector::iterator it = replaced_hydrogens.begin(); it != replaced_hydrogens.end(); it++)
    {
        PdbPreprocessorReplacedHydrogen* replaced_hydrogen = (*it);
        PdbFileSpace::PdbAtomCard* pdb_atom =
                new PdbFileSpace::PdbAtomCard(replaced_hydrogen->GetResidueChainId(),replaced_hydrogen->GetAtomName(),
                            replaced_hydrogen->GetResidueName(), replaced_hydrogen->GetResidueSequenceNumber(), replaced_hydrogen->GetResidueInsertionCode(), replaced_hydrogen->GetResidueAlternateLocation());
//        this->GetPdbFile().DeleteAtomWithTheGivenModelNumber(pdb_atom, model_number);
        to_be_deleted_atoms_.push_back(pdb_atom);
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(model_number);
}

bool PdbPreprocessor::ExtractAminoAcidChains()
{
    try
    {
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResiduesFromAtomSection();

        PdbPreprocessorChainIdResidueMap all_chain_map_residue;
        PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
        PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
        {
            PdbFileSpace::PdbResidue* residue = *it;
            char chain_id = residue->GetResidueChainId();

            std::stringstream ss;
            ss << chain_id;
            all_chain_map_residue[ss.str()].push_back(residue);
        }

        // Slower version
//        std::vector<std::string> amino_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
        // Advanced version
        gmml::ResidueNameMap amino_lib_residue_names_map = this->GetParameters().GetAllResidueNameMap();
        amino_lib_residue_names_map["HIS"] = "HIS";

        for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
        {
            int internal_amino_acid_chain_counter = 0;
            PdbPreprocessorChainIdResidueMap chain_map_residue = PdbPreprocessorChainIdResidueMap();
            std::string chain_id = (*it).first;
            PdbFileSpace::PdbFile::PdbResidueVector residues = (*it).second;

            for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
            {
                PdbFileSpace::PdbResidue* residue = (*it1);
                if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
                {
                    std::stringstream ss;
                    ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    std::stringstream sss(chain_id);
                    char c_id;
                    sss >> c_id;
                    chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                    chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
                }
                else
                {
                    std::stringstream ss;
                    ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    internal_amino_acid_chain_counter++;
                }
            }
            if(chain_map_residue.size() > 2)
            {
//                std::cout << "There is an undefined protein in the middle of the chain" << std::endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is an undefined protein in the middle of the chain" );
//                std::cout << "Pdb file is not processible at this time" << std::endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "Pdb file is not processible at this time" );


                return false;
            }
        }
        for(PdbPreprocessorChainIdSequenceNumbersMap::iterator it = chain_map_sequence_number.begin(); it != chain_map_sequence_number.end(); it++)
        {
            char chain_id = (*it).first;
            std::vector<int> sequence_numbers = (*it).second;
            std::vector<char> insertion_codes = chain_map_insertion_code[chain_id];

            std::vector<int>::iterator starting_sequence_number_iterator = min_element(sequence_numbers.begin(), sequence_numbers.end());
            std::vector<int>::iterator ending_sequence_number_iterator = max_element(sequence_numbers.begin(), sequence_numbers.end());
            int starting_index = distance(sequence_numbers.begin(), starting_sequence_number_iterator);
            int ending_index = distance(sequence_numbers.begin(), ending_sequence_number_iterator);

            PdbPreprocessorChainTermination* chain = new PdbPreprocessorChainTermination(chain_id, *starting_sequence_number_iterator, *ending_sequence_number_iterator,
                                                                                         insertion_codes.at(starting_index), insertion_codes.at(ending_index) );

            this->AddChainTermination(chain);
        }
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {
        return false;
    }
}

void PdbPreprocessor::UpdateAminoAcidChainsWithTheGivenModelNumber(PdbPreprocessorChainTerminationVector chain_terminations, int model_number)
{
    // Before non-amino-acid residue
    for(PdbPreprocessor::PdbPreprocessorChainTerminationVector::iterator it1 = chain_terminations.begin(); it1 != chain_terminations.end(); it1++)
    {
        PdbPreprocessorChainTermination* chain = (*it1);
        this->GetPdbFile().SplitAtomCardOfModelCard(chain->GetResidueChainId(), chain->GetEndingResidueSequenceNumber() + 1);
        std::stringstream ss;
        ss <<  chain->GetResidueChainId() << "_" << chain->GetEndingResidueSequenceNumber();
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Split atom card after: " + ss.str());
    }
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Putting TER card after non-amino acid residues aka chain_terminations: Done" );

    // Get all TER card positions and split
    std::vector<std::pair<char, int> > ter_card_positions = this->GetPdbFile().GetAllTerCardPositions(this->GetParameters().GetGlycamResidueNames());

    for(std::vector<std::pair<char, int> >::iterator it1 = ter_card_positions.begin(); it1 != ter_card_positions.end(); it1++)
    {
        std::pair<char, int> ter_position = *it1;
        char chain_id = ter_position.first;
        int sequence_number = ter_position.second;
        this->GetPdbFile().SplitAtomCardOfModelCard(chain_id, sequence_number);
    }
//    std::cout << "Putting TER card after residues with no tail or with at least 2 tails: Done" << std::endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Putting TER card after residues with no tail or with at least 2 tails: Done" );

    // Add Terminals
    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
    for(PdbPreprocessor::PdbPreprocessorChainTerminationVector::iterator it1 = chain_terminations.begin(); it1 != chain_terminations.end(); it1++)
    {
        PdbPreprocessorChainTermination* chain = (*it1);
        // Zwitterionic in n terminal
        if(chain->GetStringFormatOfSelectedNTermination().find("+") != std::string::npos || chain->GetStringFormatOfSelectedNTermination().find("-") != std::string::npos)
        {
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Did nothing for N termination, which selected is " + chain->GetStringFormatOfSelectedNTermination() );
            // First find out name of residue so can look up lib file
            // See if any atoms are missing, or any are extra.
            // Grow missing ones? That doesn't look easy actually, I'm not even sure there is a single point of truth for the atom cards.
        }
        else // this adds an extra residue, the above if should add hydrogens to the residue to make it charged.
        {   // Add n terminal residue at the beginning of the chain
            gmml::PossibleNChainTermination n_termination = chain->GetSelectedNTermination();
            std::string string_n_termination = chain->GetStringFormatOfNTermination(n_termination);
            LibraryFileSpace::LibraryFileResidue* lib_file_residue = this->GetParameters().FindLibResidue(string_n_termination);
            if(lib_file_residue != NULL)
            {
                LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                PdbFileSpace::PdbAtomSection* pdb_atom_card = new PdbFileSpace::PdbAtomSection();
                pdb_atom_card->SetRecordName("ATOM");
                PdbFileSpace::PdbAtomSection::PdbAtomMap atom_map;
                PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector();
                int serial_number = 0;
                for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                {
                    LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                    PdbFileSpace::PdbAtomCard* pdb_atom = new PdbFileSpace::PdbAtomCard(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, chain->GetResidueChainId(),
                            chain->GetStartingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), gmml::dNotSet, gmml::dNotSet, " ", "");
                    atom_map[serial_number] = pdb_atom;
                    atom_vector.push_back(pdb_atom);
                    serial_number++;
                }
                pdb_atom_card->SetAtomCards(atom_map);
                pdb_atom_card->SetOrderedAtomCards(atom_vector);
                gmml::log(__LINE__, __FILE__,  gmml::INF, "CIs this code ever triggered? ");
                this->GetPdbFile().InsertResidueBeforeWithTheGivenModelNumber(pdb_atom_card, model_number);
            }
        }
        // Zwitterionic in c terminal
        if(chain->GetStringFormatOfSelectedCTermination().find("+") != std::string::npos || chain->GetStringFormatOfSelectedCTermination().find("-") != std::string::npos)
        {
            // OG do nothing. TLEAP needs to add hydrogens if they are missing, as well as OXT atoms.
            // Figuring out the coords for missing atoms is too hard. Tleap is good at it.
//            OLIVER HERE
//            PdbFileSpace::PdbFile::PdbResidueVector pdbResidues = this->GetPdbFile().GetAllResiduesFromAtomSection();
//            PdbFileSpace::PdbResidue* firstRes = pdbResidues.front();
//            PdbFileSpace::PdbResidue* finalRes = pdbResidues.back();
//            std::stringstream ss;
//            ss << "firstRes:" << firstRes->GetResidueName() << firstRes->GetResidueSequenceNumber() << "." << firstRes->GetResidueChainId();
//            ss << "finalRes:" << finalRes->GetResidueName() << finalRes->GetResidueSequenceNumber() << "." << finalRes->GetResidueChainId();
//            gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str() );
//            LibraryFileSpace::LibraryFileResidue* lib_file_residue = this->GetParameters().FindLibResidue("C" + finalRes->GetResidueName());
//            //PdbAtomCardVector GetAllAtomsOfResidue(PdbResidue* residue);
//            std::vector<PdbFileSpace::PdbAtomCard*> finalResidueAtoms = this->GetPdbFile().GetAllAtomsOfResidue(finalRes);
//            // Go through lib file atoms, if one exists that isn't in residue, add it.
//            for (auto &libFileAtom : lib_file_residue->GetAtomsVector())
//            {
//                bool found = false;
//                for (auto &pdbAtomCard : finalResidueAtoms)
//                {
//                    if (pdbAtomCard->GetAtomName() == libFileAtom->GetName())
//                    {
//                        found = true;
//                    }
//                }
//                if (!found)
//                {
//                    gmml::log(__LINE__, __FILE__,  gmml::INF, "Creating a new atom called " + libFileAtom->GetName());
//                                                            int serial_number = 900000;
//                                                            GeometryTopology::Coordinate coords(0.0, 0.0, 0.0);
            //GeometryTopology::Coordinate cCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(cCoordProtein, caCoordProtein, nCoordProtein, 120.0, -130.0, 1.4);
            //                    PdbAtomCard(
            //                    int atom_serial_number, std::string atom_name, char atom_alternate_location, std::string residue_name, char chain_id,
            //                                        int residue_sequence_number, char insertion_code, GeometryTopology::Coordinate coordinate, double occupancy, double tempreture_factor,
            //                                        std::string element_symbol, std::string charge, std::vector<PdbAtomCard*> alternate_atom_locations ={});
            //                    PdbFileSpace::PdbAtomCard new_atom(
            //                            serial_number, //  int atom_serial_number
            //                            libFileAtom->GetName(), // std::string atom_name
            //                            firstRes->GetResidueAlternateLocation(), //char atom_alternate_location
            //                            firstRes->GetResidueName(), // std::string residue_name,
            //                            firstRes->GetResidueChainId(), // char chain_id
            //                            firstRes->GetResidueSequenceNumber(), // int residue_sequence_number,
            //                            firstRes->GetResidueInsertionCode(), // char insertion_code,
            //                            coords, // GeometryTopology::Coordinate coordinate
            //                            1.00, //  double occupancy
            //                            0.00, // double tempreture_factor
            //                            std::string(), // std::string element_symbol
            //                            std::string(), //  std::string charge
            //                            std::vector<PdbFileSpace::PdbAtomCard*> {}
            //                    );
//                }
//            }
//
            // Go through residue atoms, if one exists that isn't in the lib, delete it or cry idk.
        }
        else
        {
            // Add c terminal at the end of the chain
            gmml::PossibleCChainTermination c_termination = chain->GetSelectedCTermination();
            std::string string_c_termination = chain->GetStringFormatOfCTermination(c_termination);
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Adding c termination: " + string_c_termination );
            LibraryFileSpace::LibraryFileResidue* lib_file_residue = this->GetParameters().FindLibResidue(string_c_termination);
            if(lib_file_residue != NULL)
            {
                lib_atoms = lib_file_residue->GetAtoms();
                gmml::log(__LINE__, __FILE__,  gmml::INF, "Found lib file: " + lib_file_residue->GetName() );
                PdbFileSpace::PdbAtomSection* pdb_atom_card = new PdbFileSpace::PdbAtomSection();
                pdb_atom_card->SetRecordName("ATOM");
                PdbFileSpace::PdbAtomSection::PdbAtomMap atom_map;
                PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector();
                int serial_number = 0;
                for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                {
                    LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                    PdbFileSpace::PdbAtomCard* pdb_atom = new PdbFileSpace::PdbAtomCard(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, chain->GetResidueChainId(),
                            chain->GetEndingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), gmml::dNotSet, gmml::dNotSet, " ", "");
                    atom_map[serial_number] = pdb_atom;
                    atom_vector.push_back(pdb_atom);
                    serial_number++;
                }
                pdb_atom_card->SetAtomCards(atom_map);
                pdb_atom_card->SetOrderedAtomCards(atom_vector);
                gmml::log(__LINE__, __FILE__,  gmml::INF, "AIs this code ever triggered? ");
                this->GetPdbFile().InsertResidueAfterWithTheGivenModelNumber(pdb_atom_card, model_number);
            }
        }
    }
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Add N/C terminals: Done" );
    return;
}

bool PdbPreprocessor::ExtractGapsInAminoAcidChains()
{
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResiduesFromAtomSection();

    PdbPreprocessorChainIdResidueMap all_chain_map_residue;
    PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
    PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* residue = *it;
        char chain_id = residue->GetResidueChainId();

        std::stringstream ss;
        ss << chain_id;
        all_chain_map_residue[ss.str()].push_back(residue);
    }

    // Slower version
//        std::vector<std::string> amino_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
    // Advanced version
    gmml::ResidueNameMap amino_lib_residue_names_map = this->GetParameters().GetAllResidueNameMap(); // not exactly the same, returns more than just the amino.
    amino_lib_residue_names_map["HIS"] = "HIS";
    PdbPreprocessorChainIdResidueMap all_chain_map_amino_acid_residue = PdbPreprocessorChainIdResidueMap();

    PdbFileSpace::PdbFile::PdbResidueVector residues;
    PdbPreprocessorChainIdResidueMap chain_map_residue;
    for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
    {
        int internal_amino_acid_chain_counter = 0;
        chain_map_residue = PdbPreprocessorChainIdResidueMap();
        std::string chain_id = (*it).first;
        residues = (*it).second;

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            PdbFileSpace::PdbResidue* residue = (*it1);
            if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
            {
                std::stringstream ss;
                ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                chain_map_residue[ss.str()].push_back(residue);
                all_chain_map_amino_acid_residue[chain_id].push_back(residue);
                std::stringstream sss(chain_id);
                char c_id;
                sss >> c_id;
                chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
            }
            else
            {
                std::stringstream ss;
                ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                chain_map_residue[ss.str()].push_back(residue);
                internal_amino_acid_chain_counter++;
            }
        }
        if(chain_map_residue.size() > 2)
        {
//            std::cout << "There is an undefined protein in the middle of the chain" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is an undefined protein in the middle of the chain" );
//            std::cout << "Pdb file is not processible at this time" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Pdb file is not processible at this time" );

            return false;
        }
    }


    PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atom_map;
    residue_atom_map = this->GetPdbFile().GetAllAtomsOfResidues();
    for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_amino_acid_residue.begin(); it != all_chain_map_amino_acid_residue.end(); it++)
    {
        std::string chain_id = (*it).first;
        PdbFileSpace::PdbFile::PdbResidueVector residues = (*it).second;
        char c_id;
        std::stringstream ss(chain_id);
        ss >> c_id;
        std::vector<int> sequence_numbers = chain_map_sequence_number[c_id];
        std::vector<char> insertion_codes = chain_map_insertion_code[c_id];
        std::vector<int>::iterator starting_sequence_number_iterator = min_element(sequence_numbers.begin(), sequence_numbers.end());
        std::vector<int>::iterator ending_sequence_number_iterator = max_element(sequence_numbers.begin(), sequence_numbers.end());

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            unsigned int dist = distance(residues.begin(), it1);
            if(dist != residues.size() - 1)
            {
                PdbFileSpace::PdbResidue* residue = (*it1);
                PdbFileSpace::PdbResidue* next_residue = *(it1 + 1);
                int i = distance(residues.begin(), it1);
                int j = distance(residues.begin(), it1 + 1);
                PdbFileSpace::PdbAtomCard* c_atom_of_residue = this->GetPdbFile().GetAtomOfResidueByName(residue, "C", residue_atom_map);
                PdbFileSpace::PdbAtomCard* n_atom_of_next_residue = this->GetPdbFile().GetAtomOfResidueByName(next_residue, "N", residue_atom_map);
                double distance = 0.0;
                if(c_atom_of_residue != NULL && n_atom_of_next_residue != NULL)
                {
                    GeometryTopology::Coordinate c_atom_coordinate = c_atom_of_residue->GetAtomOrthogonalCoordinate();
                    GeometryTopology::Coordinate n_atom_coordinate = n_atom_of_next_residue->GetAtomOrthogonalCoordinate();
                    distance = c_atom_coordinate.Distance(n_atom_coordinate);
                }
                if(distance > gmml::maxCutOff + 1.0)
                {
                    PdbPreprocessorMissingResidue* missing_residues = new PdbPreprocessorMissingResidue(c_id, *starting_sequence_number_iterator,
                                                                                                        *ending_sequence_number_iterator, sequence_numbers.at(i),
                                                                                                        sequence_numbers.at(j), insertion_codes.at(i), insertion_codes.at(j));
                    this->AddMissingResidue(missing_residues);
                }
            }
        }
    }

    return true;
}

void PdbPreprocessor::UpdateGapsInAminoAcidChainsWithTheGivenModelNumber(PdbPreprocessorMissingResidueVector gaps, int model_number)
{
    for(PdbPreprocessor::PdbPreprocessorMissingResidueVector::iterator it1 = gaps.begin(); it1 != gaps.end(); it1++)
    {
        PdbPreprocessorMissingResidue* gap = (*it1);
        this->GetPdbFile().SplitAtomCardOfModelCardWithTheGivenModelNumber(gap->GetResidueChainId(), gap->GetResidueAfterGap(), model_number);
    }
    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
    for(PdbPreprocessor::PdbPreprocessorMissingResidueVector::iterator it1 = gaps.begin(); it1 != gaps.end(); it1++)
    {
        PdbPreprocessorMissingResidue* gap = (*it1);
        // Zwitterionic in n terminal
        if(gap->GetStringFormatOfSelectedNTermination().find("+") != std::string::npos || gap->GetStringFormatOfSelectedNTermination().find("-") != std::string::npos)
        {
            // Zwitterionic in c terminal
            if(gap->GetStringFormatOfSelectedCTermination().find("+") != std::string::npos || gap->GetStringFormatOfSelectedCTermination().find("-") != std::string::npos)
            {

            }
            else
            {
                // Add c terminal at the end of the chain
                gmml::PossibleCChainTermination c_termination = gap->GetSelectedCTermination();
                std::string string_c_termination = gap->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = this->GetParameters().FindLibResidue(string_c_termination);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomSection* pdb_atom_card = new PdbFileSpace::PdbAtomSection();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomSection::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtomCard* pdb_atom = new PdbFileSpace::PdbAtomCard(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueBeforeGap(), ' ', lib_file_atom->GetCoordinate(), gmml::dNotSet, gmml::dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtomCards(atom_map);
                    pdb_atom_card->SetOrderedAtomCards(atom_vector);
                    this->GetPdbFile().InsertResidueAfterWithTheGivenModelNumber(pdb_atom_card, model_number);
                }
            }
        }
        else
        {
            // Zwitterionic in c terminal
            if(gap->GetStringFormatOfSelectedCTermination().find("+") != std::string::npos || gap->GetStringFormatOfSelectedCTermination().find("-") != std::string::npos)
            {
                // Add n terminal residue at the beginning of the chain
                gmml::PossibleNChainTermination n_termination = gap->GetSelectedNTermination();
                std::string string_n_termination = gap->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = this->GetParameters().FindLibResidue(string_n_termination);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomSection* pdb_atom_card = new PdbFileSpace::PdbAtomSection();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomSection::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtomCard* pdb_atom = new PdbFileSpace::PdbAtomCard(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueAfterGap(), ' ', lib_file_atom->GetCoordinate(), gmml::dNotSet, gmml::dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtomCards(atom_map);
                    pdb_atom_card->SetOrderedAtomCards(atom_vector);
                    this->GetPdbFile().InsertResidueBeforeWithTheGivenModelNumber(pdb_atom_card, model_number);
                    if(it1 != --gaps.end())
                    {
                        PdbPreprocessorMissingResidue* next_gap = (*(++it1));
                        if(gap->GetResidueChainId() == next_gap->GetResidueChainId())
                        {
                            (*it1)->SetResidueBeforeGap((*it1)->GetResidueBeforeGap()+1);
                        }
                        it1--;
                    }
                }
            }
            else
            {
                // Add n terminal residue at the beginning of the chain
                gmml::PossibleNChainTermination n_termination = gap->GetSelectedNTermination();
                std::string string_n_termination = gap->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = this->GetParameters().FindLibResidue(string_n_termination);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomSection* pdb_atom_card = new PdbFileSpace::PdbAtomSection();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomSection::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector = PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtomCard* pdb_atom = new PdbFileSpace::PdbAtomCard(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueAfterGap(), ' ', lib_file_atom->GetCoordinate(), gmml::dNotSet, gmml::dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtomCards(atom_map);
                    pdb_atom_card->SetOrderedAtomCards(atom_vector);
                    this->GetPdbFile().InsertResidueBeforeWithTheGivenModelNumber(pdb_atom_card, model_number);
                }
                // Add c terminal residue at the end of the chain
                gmml::PossibleCChainTermination c_termination = gap->GetSelectedCTermination();
                std::string string_c_termination = gap->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue_from_c_termination = this->GetParameters().FindLibResidue(string_c_termination);
                if(lib_file_residue_from_c_termination != nullptr)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms_from_c_termination = lib_file_residue_from_c_termination->GetAtoms();
                    PdbFileSpace::PdbAtomSection* pdb_atom_card_for_c_termination = new PdbFileSpace::PdbAtomSection();
                    pdb_atom_card_for_c_termination->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomSection::PdbAtomMap atom_map_for_c_termination;
                    PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector atom_vector_for_c_termination = PdbFileSpace::PdbAtomSection::PdbAtomCardOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it3 = lib_atoms_from_c_termination.begin(); it3 != lib_atoms_from_c_termination.end(); it3++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it3).second;
                        PdbFileSpace::PdbAtomCard* pdb_atom = new PdbFileSpace::PdbAtomCard(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueBeforeGap(), ' ', lib_file_atom->GetCoordinate(), gmml::dNotSet, gmml::dNotSet, " ", "");
                        atom_map_for_c_termination[serial_number] = pdb_atom;
                        atom_vector_for_c_termination.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card_for_c_termination->SetAtomCards(atom_map_for_c_termination);
                    pdb_atom_card_for_c_termination->SetOrderedAtomCards(atom_vector_for_c_termination);
                    this->GetPdbFile().InsertResidueAfterWithTheGivenModelNumber(pdb_atom_card_for_c_termination, model_number);
                    if(it1 != --gaps.end())
                    {
                        PdbPreprocessorMissingResidue* next_gap = (*(++it1));
                        if(gap->GetResidueChainId() == next_gap->GetResidueChainId())
                        {
                            (*it1)->SetResidueBeforeGap((*it1)->GetResidueBeforeGap()+1);
                        }
                        it1--;
                    }
                }
            }
        }
    }
}


bool PdbPreprocessor::ExtractAlternateResidue()
{
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResidues();
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* target_residue = *it;
        std::string target_residue_name = target_residue->GetResidueName();
        char target_chain_id = target_residue->GetResidueChainId();
        int target_sequence_number = target_residue->GetResidueSequenceNumber();
        char target_insertion_code = target_residue->GetResidueInsertionCode();
        char target_alternate_location = target_residue->GetResidueAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_chain_id << "_" << target_sequence_number << "_" << target_insertion_code;
        std::string target_key = ss.str();
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 =it+1 ; it1 != pdb_residues.end(); it1++)
        {
            PdbFileSpace::PdbResidue* residue = *it1;
            std::string residue_name = residue->GetResidueName();
            char chain_id = residue->GetResidueChainId();
            int sequence_number = residue->GetResidueSequenceNumber();
            char insertion_code = residue->GetResidueInsertionCode();
            char alternate_location = residue->GetResidueAlternateLocation();
            std::stringstream sss;
            sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code;
            std::string key = sss.str();

            if(key.compare(target_key) == 0)
            {
                if(target_alternate_location != alternate_location)
                {
                    if (alternate_residue_map_.empty() || distance(alternate_residue_map_.begin() ,alternate_residue_map_.find(key)) < 0 ||
                            distance(alternate_residue_map_.begin() ,alternate_residue_map_.find(key)) >= (int)alternate_residue_map_.size())
                    {
                        std::vector<bool> selected = std::vector<bool>();
                        selected.push_back(true);
                        selected.push_back(true); // OG edit, changed to true so that A coords get written out
                        std::vector<char> alternate_locations = std::vector<char>();
                        alternate_locations.push_back(target_alternate_location);
                        alternate_locations.push_back(alternate_location);
                        alternate_residue_map_[target_key] = new PdbPreprocessorAlternateResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_locations, selected);
                    }
                    else
                    {
                        PdbPreprocessorAlternateResidue* alternate_residue = alternate_residue_map_[target_key];
                        std::vector<char> alternate_locations = alternate_residue->GetResidueAlternateLocation();
                        if(distance(alternate_locations.begin(), find(alternate_locations.begin(),alternate_locations.end(),alternate_location)) < 0 ||
                                distance(alternate_locations.begin(), find(alternate_locations.begin(),alternate_locations.end(),alternate_location)) >= (int)alternate_locations.size())
                        {

                            alternate_locations.push_back(alternate_location);
                            alternate_residue->SetResidueAlternateLocation(alternate_locations);
                            std::vector<bool> selected_alternate_locations = alternate_residue->GetSelectedAlternateLocation();
                            selected_alternate_locations.push_back(false);
                            alternate_residue->SetSelectedAlternateLocation(selected_alternate_locations);
                        }
                    }
                }
            }
        }
    }
    return true;
}

void PdbPreprocessor::RemoveUnselectedAlternateResiduesWithTheGivenModelNumber(PdbPreprocessorAlternateResidueMap alternate_residue_map/*, int model_number*/)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorAlternateResidueMap::iterator it = alternate_residue_map.begin(); it != alternate_residue_map.end(); it++)
    {
        PdbPreprocessorAlternateResidue* alternate_residue = (*it).second;
        std::vector<bool> selected_alternate_locations = alternate_residue->GetSelectedAlternateLocation();
        for(std::vector<bool>::iterator it1 = selected_alternate_locations.begin(); it1 != selected_alternate_locations.end(); it1++)
        {
            bool selected_alternate_location = (*it1);
            char alternate_location = alternate_residue->GetResidueAlternateLocation().at(distance(selected_alternate_locations.begin(), it1));
            if(!selected_alternate_location)
            {
                PdbFileSpace::PdbResidue* pdb_residue = new PdbFileSpace::PdbResidue(alternate_residue->GetResidueName() ,alternate_residue->GetResidueChainId(), alternate_residue->GetResidueSequenceNumber(),
                                                         alternate_residue->GetResidueInsertionCode(), alternate_location);
//                this->GetPdbFile().DeleteResidueWithTheGivenModelNumber(pdb_residue, model_number);
                to_be_deleted_residues_.push_back(pdb_residue);
            }
        }
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber();
}

void PdbPreprocessor::DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(int model_number)
{
    if(to_be_deleted_atoms_.size() != 0)
        this->GetPdbFile().DeleteAtomsWithTheGivenModelNumber(to_be_deleted_atoms_, model_number);

    if(to_be_deleted_residues_.size() != 0)
        this->GetPdbFile().DeleteResiduesWithTheGivenModelNumber(to_be_deleted_residues_, model_number);
}

bool PdbPreprocessor::ExtractResidueInfo()
{
    LibraryFileSpace::LibraryFile::ResidueMap lib_residues = this->GetParameters().GetLibraryResidueMap();
    PrepFileSpace::PrepFile::ResidueMap prep_residues = this->GetParameters().GetPrepResidueMap();
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = this->GetPdbFile().GetAllResidues();

    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
    PrepFileSpace::PrepFileAtomVector prep_atoms;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        double residue_charge = 0.0;
        PdbFileSpace::PdbResidue* residue = *it;
        if(lib_residues.find(residue->GetResidueName()) != lib_residues.end())
        {
            LibraryFileSpace::LibraryFileResidue* lib_residue = lib_residues[residue->GetResidueName()];
            lib_atoms = lib_residue->GetAtoms();
            for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it1 = lib_atoms.begin(); it1 != lib_atoms.end(); it1++)
            {
                LibraryFileSpace::LibraryFileAtom* lib_atom = (*it1).second;
                residue_charge += lib_atom->GetCharge();
            }
        }
        else if (prep_residues.find(residue->GetResidueName()) != prep_residues.end())
        {
            PrepFileSpace::PrepFileResidue* prep_residue = prep_residues[residue->GetResidueName()];
            prep_atoms = prep_residue->GetAtoms();
            for(PrepFileSpace::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
            {
                PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
                residue_charge += prep_atom->GetCharge();
            }
        }
        std::stringstream ss;
        ss << residue->GetResidueName() << "_" << residue->GetResidueChainId() << "_" << residue->GetResidueSequenceNumber() << "_" << residue->GetResidueInsertionCode()
           << "_" << residue->GetResidueAlternateLocation();
        PdbPreprocessorResidueInfo* residue_info = new PdbPreprocessorResidueInfo(residue->GetResidueName(), residue->GetResidueChainId(), residue->GetResidueSequenceNumber(),
                                                                                  residue->GetResidueInsertionCode(), residue->GetResidueAlternateLocation(), residue_charge);
        residue_info_map_[ss.str()] = residue_info;
    }
    return true;
}

double PdbPreprocessor::CalculateModelCharge()
{
    double modelCharge = 0.0;
    for(auto &pdbRes : this->GetPdbFile().GetAllResidues())
    {
        modelCharge += this->GetParameters().GetChargeForResidue(pdbRes->GetResidueName());
    }
    return modelCharge;
}

//void PdbPreprocessor::Preprocess(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files_path, std::vector<std::string> glycam_lib_files_path, std::vector<std::string> other_lib_files_path, std::vector<std::string> prep_files_path)
void PdbPreprocessor::Preprocess()
{
    try
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, " Start preprocessing ..." );
        if(this->ExtractHISResidues())
            gmml::log(__LINE__, __FILE__,  gmml::INF, " HIS residues extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF, " HIS residues extraction: failed");
        if(this->ExtractCYSResidues())
            gmml::log(__LINE__, __FILE__,  gmml::INF, " CYS residues extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF, " CYS residues extraction: failed");
        if(this->ExtractAlternateResidue())
            gmml::log(__LINE__, __FILE__,  gmml::INF, " Alternate residues extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF, " Alternate residues extraction: failed");
        if(this->ExtractUnrecognizedResidues())
            gmml::log(__LINE__, __FILE__,  gmml::INF, " Unrecognized residues extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF, " Unrecognized residues extraction: failed");
        if(this->ExtractUnknownHeavyAtoms())
            gmml::log(__LINE__, __FILE__,  gmml::INF, " Unknown heavy atoms extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF,  " Unknown heavy atoms extraction: failed");
        if(this->ExtractRemovedHydrogens())
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Removed hydrogens extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Removed hydrogens extraction: failed");
        if(this->ExtractAminoAcidChains())
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Amino acid chains extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Amino acid chains extraction: failed");
        if(this->ExtractGapsInAminoAcidChains())
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Gaps in amino acid chains extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Gaps in amino acid chains extraction: failed");
        if(this->ExtractResidueInfo())
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Residue info extraction: done");
        else
            gmml::log(__LINE__, __FILE__,  gmml::INF, "Residue info extraction: failed");
        std::stringstream model_charge;
        model_charge << "Model charge is " << this->CalculateModelCharge() ;
        gmml::log(__LINE__, __FILE__,  gmml::INF, model_charge.str() );
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Model charge calculation: done");
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Preprocessing done");
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Caught an exception during Preprocessing: " + std::string(ex.what()));
    }
}

void PdbPreprocessor::ApplyPreprocessingWithTheGivenModelNumber(int model_number)
{
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Start to apply changes ...");
    UpdateHISMappingWithTheGivenNumber(this->GetHistidineMappings(), model_number);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "HIS residues update: done");
    UpdateCYSResiduesWithTheGivenModelNumber(this->GetDisulfideBonds());
    gmml::log(__LINE__, __FILE__,  gmml::INF, "CYS residues update: done");
    RemoveUnselectedAlternateResiduesWithTheGivenModelNumber(this->GetAlternateResidueMap()/*, model_number*/);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Unselected alternate residues removed: done");
    RemoveUnrecognizedResiduesWithTheGivenModelNumber(this->GetUnrecognizedResidues(), model_number);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Remove unrecognized residues: done");
    RemoveResiduesOfUnknownHeavyAtomsWithTheGivenModelNumber(this->GetUnrecognizedHeavyAtoms(), model_number);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Unknown heavy atoms removed: done");
    RemoveRemovedHydrogensWithTheGivenModelNumber(this->GetReplacedHydrogens(), model_number);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Removed hydrogens removed: done");
    UpdateAminoAcidChainsWithTheGivenModelNumber(this->GetChainTerminations(), model_number);
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Amino acid chains update: done");
    UpdateGapsInAminoAcidChainsWithTheGivenModelNumber(this->GetMissingResidues(), model_number); // OG Mar 2017
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Fixed the gaps.");
    std::stringstream model_charge;
    model_charge << "Final model charge is " << this->CalculateModelCharge() ;
    gmml::log(__LINE__, __FILE__,  gmml::INF, model_charge.str() );
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Applying changes completed.");
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessor::Print(std::ostream &out)
{
    out << "================================== Disulfide Bonds =====================================" << std::endl;
    for(auto &disulfide_bond : this->GetDisulfideBonds())
    {
        disulfide_bond->Print(out);
    }
    out << "======================================= Chains =========================================" << std::endl;
    for(auto &chain : this->GetChainTerminations())
    {
        chain->Print(out);
    }
    out << "===================================== HIS Residues =====================================" << std::endl;
    for(auto &hisResidue : this->GetHistidineMappings())
    {
        hisResidue->Print(out);
    }
    out << "======================================== Gaps ==========================================" << std::endl;
    for(auto &gap : this->GetMissingResidues())
    {
        gap->Print(out);
    }
    out << "============================== Unrecognized Residues ===================================" << std::endl;
    for(auto &unrecognizedResidue : this->GetUnrecognizedResidues())
    {
        unrecognizedResidue->Print(out);
    }
    out << "=============================== Unknown Heavy Atoms ====================================" << std::endl;
    for(auto &unknownHeavyAtom : this->GetUnrecognizedHeavyAtoms())
    {
        unknownHeavyAtom->Print(out);
    }
    out << "================================ Removed Hydrogens =====================================" << std::endl;
    for(auto &removedHydrogen : this->GetReplacedHydrogens())
    {
        removedHydrogen->Print(out);
    }
    out << "================================ Alternate Residues ====================================" << std::endl;
    for(auto &alternateResidue : this->GetAlternateResidueMap())
    {
        alternateResidue.second->Print(out);
    }
    out << "================================ Residues Info ====================================" << std::endl;
    for(auto &residueInfo : this->GetResidueInfoMap())
    {
        residueInfo.second->Print(out);
    }
}
