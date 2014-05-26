#include <cctype>

#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessor.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp"
#include "../../../includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/common.hpp";
#include "../../../includes/utils.hpp"
using namespace std;
using namespace PdbPreprocessorSpace;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbPreprocessor::PdbPreprocessor() {}

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
PdbPreprocessor::PdbPreprocessorRecognizedResidueVector PdbPreprocessor::GetRecognizedResidues(){
    return recognized_residues_;
}
PdbPreprocessor::PdbPreprocessorUnrecognizedHeavyAtomVector PdbPreprocessor::GetUnrecognizedHeavyAtoms(){
    return unrecognized_heavy_atoms_;
}
PdbPreprocessor::PdbPreprocessorReplacedHydrogenVector PdbPreprocessor::GetReplacedHydrogens(){
    return replaced_hydrogens_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbPreprocessor::SetDisulfideBonds(PdbPreprocessor::PdbPreprocessorDisulfideBondVector disulfide_bonds){
    disulfide_bonds_.clear();
    for(PdbPreprocessorDisulfideBondVector::iterator it = disulfide_bonds.begin(); it != disulfide_bonds.end(); it++)
    {
        disulfide_bonds_.push_back(*it);
    }
}
void PdbPreprocessor::AddDisulfideBond(PdbPreprocessorDisulfideBond *disulfide_bond)
{
    disulfide_bonds_.push_back(disulfide_bond);
}
void PdbPreprocessor::SetChainTerminations(PdbPreprocessor::PdbPreprocessorChainTerminationVector chain_terminations){
    chain_terminations_.clear();
    for(PdbPreprocessorChainTerminationVector::iterator it = chain_terminations.begin(); it != chain_terminations.end(); it++)
    {
        chain_terminations_.push_back(*it);
    }
}
void PdbPreprocessor::AddChainTermination(PdbPreprocessorChainTermination *chain_termination)
{
    chain_terminations_.push_back(chain_termination);
}
void PdbPreprocessor::SetHistidineMappings(PdbPreprocessor::PdbPreprocessorHistidineMappingVector histidine_mappings){
    histidine_mappings_.clear();
    for(PdbPreprocessorHistidineMappingVector::iterator it = histidine_mappings.begin(); it != histidine_mappings.end(); it++)
    {
        histidine_mappings_.push_back(*it);
    }
}
void PdbPreprocessor::AddHistidineMapping(PdbPreprocessorHistidineMapping *histidine_mapping)
{
    histidine_mappings_.push_back(histidine_mapping);
}
void PdbPreprocessor::SetMissingResidues(PdbPreprocessor::PdbPreprocessorMissingResidueVector missing_residues){
    missing_residues_.clear();
    for(PdbPreprocessorMissingResidueVector::iterator it = missing_residues.begin(); it != missing_residues.end(); it++)
    {
        missing_residues_.push_back(*it);
    }
}
void PdbPreprocessor::AddMissingResidue(PdbPreprocessorMissingResidue *missing_residue)
{
    missing_residues_.push_back(missing_residue);
}
void PdbPreprocessor::SetUnrecognizedResidues(PdbPreprocessor::PdbPreprocessorUnrecognizedResidueVector unrecognized_residues){
    unrecognized_residues_.clear();
    for(PdbPreprocessorUnrecognizedResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
    {
        unrecognized_residues_.push_back(*it);
    }
}
void PdbPreprocessor::AddUnrecognizedResidue(PdbPreprocessorUnrecognizedResidue *unrecognized_residue)
{
    unrecognized_residues_.push_back(unrecognized_residue);
}
void PdbPreprocessor::SetRecognizedResidues(PdbPreprocessor::PdbPreprocessorRecognizedResidueVector recognized_residues){
    recognized_residues_.clear();
    for(PdbPreprocessorRecognizedResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        recognized_residues_.push_back(*it);
    }
}
void PdbPreprocessor::AddRecognizedResidue(PdbPreprocessorUnrecognizedResidue *recognized_residue)
{
    recognized_residues_.push_back(recognized_residue);
}
void PdbPreprocessor::SetUnrecognizedHeavyAtoms(PdbPreprocessor::PdbPreprocessorUnrecognizedHeavyAtomVector unrecognized_heavy_atoms){
    unrecognized_heavy_atoms_.clear();
    for(PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unrecognized_heavy_atoms.begin(); it != unrecognized_heavy_atoms.end(); it++)
    {
        unrecognized_heavy_atoms_.push_back(*it);
    }
}
void PdbPreprocessor::AddUnrecognizedHeavyAtom(PdbPreprocessorUnrecognizedHeavyAtom *unrecognized_heavy_atom)
{
    unrecognized_heavy_atoms_.push_back(unrecognized_heavy_atom);
}
void PdbPreprocessor::SetReplacedHydrogens(PdbPreprocessor::PdbPreprocessorReplacedHydrogenVector replaced_hydrogens){
    replaced_hydrogens.clear();
    for(PdbPreprocessorReplacedHydrogenVector::iterator it = replaced_hydrogens.begin(); it != replaced_hydrogens.end(); it++)
    {
        replaced_hydrogens_.push_back(*it);
    }
}
void PdbPreprocessor::AddReplacedHydrogen(PdbPreprocessorReplacedHydrogen *replaced_hydrogen)
{
    replaced_hydrogens_.push_back(replaced_hydrogen);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
vector<string> PdbPreprocessor::GetUnrecognizedResidueNames(vector<string> pdb_residue_names, vector<string> dataset_residue_names)
{
    vector<string> unrecognized_residue_names;
    bool is_recognized = false;
    for(vector<string>::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        string pdb_residue_name = *it;
        for(vector<string>::iterator it1 = dataset_residue_names.begin(); it1 != dataset_residue_names.end(); it1++)
        {
            string dataset_residue_name = *it1;
            if((pdb_residue_name).compare(dataset_residue_name) == 0)
            {
                is_recognized = true;
                break;
            }
            else
                is_recognized = false;
        }
        if(!is_recognized)
            unrecognized_residue_names.push_back(pdb_residue_name);
    }
    return unrecognized_residue_names;
}
vector<string> PdbPreprocessor::GetRecognizedResidueNames(vector<string> pdb_residue_names, vector<string> dataset_residue_names)
{
    vector<string> recognized_residue_names;
    for(vector<string>::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        string pdb_residue_name = *it;
        for(vector<string>::iterator it1 = dataset_residue_names.begin(); it1 != dataset_residue_names.end(); it1++)
        {
            string dataset_residue_name = *it1;
            if((pdb_residue_name).compare(dataset_residue_name) == 0)
            {
                recognized_residue_names.push_back(pdb_residue_name);
                break;
            }
        }
    }
    return recognized_residue_names;
}
PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetUnrecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, vector<string> unrecognized_residue_names)
{
    PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues;
    bool is_unrecognized = false;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* pdb_residue = *it;
        string pdb_residue_name = pdb_residue->GetResidueName();
        for(vector<string>::iterator it1 = unrecognized_residue_names.begin(); it1 != unrecognized_residue_names.end(); it1++)
        {
            string unrecognized_residue_name = *it1;
            if((pdb_residue_name).compare(unrecognized_residue_name) == 0)
            {
                is_unrecognized = true;
                break;
            }
            else
                is_unrecognized = false;
        }
        if(is_unrecognized)
            unrecognized_residues.push_back(pdb_residue);
    }
    return unrecognized_residues;
}
PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetRecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, vector<string> recognized_residue_names)
{
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues;
    bool is_recognized = false;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* pdb_residue = *it;
        string pdb_residue_name = pdb_residue->GetResidueName();
        for(vector<string>::iterator it1 = recognized_residue_names.begin(); it1 != recognized_residue_names.end(); it1++)
        {
            string recognized_residue_name = *it1;
            if((pdb_residue_name).compare(recognized_residue_name) == 0)
            {
                is_recognized = true;
                break;
            }
            else
                is_recognized = false;
        }
        if(is_recognized)
            recognized_residues.push_back(pdb_residue);
    }
    return recognized_residues;
}



vector<string> PdbPreprocessor::GetAllResidueNamesFromMultipleLibFiles(vector<string> lib_files)
{
    vector<string> all_residue_names;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        vector<string> residue_names = lib_file->GetAllResidueNames();
        for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
        {
            all_residue_names.push_back(*it1);
        }
    }
    return all_residue_names;
}

vector<string> PdbPreprocessor::GetAllResidueNamesFromMultiplePrepFiles(vector<string> prep_files)
{
    vector<string> all_residue_names;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        vector<string> residue_names = prep_file->GetAllResidueNames();
        for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
        {
            all_residue_names.push_back(*it1);
        }
    }
    return all_residue_names;
}

vector<string> PdbPreprocessor::GetAllResidueNamesFromDatasetFiles(vector<string> lib_files, vector<string> prep_files)
{
    vector<string> all_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
    vector<string> all_prep_residue_names = GetAllResidueNamesFromMultiplePrepFiles(prep_files);
    vector<string> all_residue_names;
    for(vector<string>::iterator it1 = all_lib_residue_names.begin(); it1 != all_lib_residue_names.end(); it1++)
    {
        all_residue_names.push_back(*it1);
    }
    for(vector<string>::iterator it1 = all_prep_residue_names.begin(); it1 != all_prep_residue_names.end(); it1++)
    {
        all_residue_names.push_back(*it1);
    }
    return all_residue_names;
}

void PdbPreprocessor::ExtractUnrecognizedResidues(string pdb_file_path, vector<string> lib_files, vector<string> prep_files)
{
    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    PdbFile* pdb_file = new PdbFile(pdb_file_path);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    vector<string> unrecognized_residue_names = GetUnrecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues = GetUnrecognizedResidues(pdb_residues, unrecognized_residue_names);

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
    {
        PdbResidue* pdb_residue = (*it);
        PdbPreprocessorUnrecognizedResidue* unrecognized_residue =
                new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(),
                                                       pdb_residue->GetResidueSequenceNumber(), pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation());
        unrecognized_residues_.push_back(unrecognized_residue);
    }
}

void PdbPreprocessor::RemoveUnrecognizedResidues(PdbFile *pdb_file, PdbPreprocessorUnrecognizedResidueVector unrecognized_residues)
{
    for(PdbPreprocessorUnrecognizedResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
    {
        PdbPreprocessorUnrecognizedResidue* unrecognized_residue = (*it);
        PdbResidue* pdb_residue =
                new PdbResidue(unrecognized_residue->GetResidueName(), unrecognized_residue->GetResidueChainId(),
                               unrecognized_residue->GetResidueSequenceNumber(), unrecognized_residue->GetResidueInsertionCode(), unrecognized_residue->GetResidueAlternateLocation());
        pdb_file->DeleteResidue(pdb_residue);
    }
}

void PdbPreprocessor::ExtractRecognizedResidues(string pdb_file_path, vector<string> lib_files, vector<string> prep_files)
{
    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    PdbFile* pdb_file = new PdbFile(pdb_file_path);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    vector<string> unrecognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, unrecognized_residue_names);

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        PdbResidue* pdb_residue = (*it);
        PdbPreprocessorUnrecognizedResidue* recognized_residue =
                new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(),
                                                       pdb_residue->GetResidueSequenceNumber(), pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation());
        recognized_residues_.push_back(recognized_residue);
    }
}

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetAllCYSResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues)
{
    PdbFileSpace::PdbFile::PdbResidueVector all_cys_residues;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* pdb_residue = (*it);
        string pdb_residue_name = pdb_residue->GetResidueName();
        if((pdb_residue_name).compare("CYS") == 0)
        {
            all_cys_residues.push_back(pdb_residue);
        }
    }
    return all_cys_residues;
}

double PdbPreprocessor::GetDistanceofCYS(PdbResidue* first_residue, PdbResidue* second_residue, PdbFile* pdb_file, PdbFile::PdbResidueAtomsMap residue_atom_map)
{
    double distance = 0.0;
    PdbAtom* first_residue_sulfur_atom = pdb_file->GetAtomOfResidueByName(first_residue, "SG", residue_atom_map);
    PdbAtom* second_residue_sulfur_atom = pdb_file->GetAtomOfResidueByName(second_residue, "SG", residue_atom_map);
    if(first_residue_sulfur_atom != NULL && second_residue_sulfur_atom != NULL)
        distance = first_residue_sulfur_atom->GetAtomOrthogonalCoordinate().Distance(second_residue_sulfur_atom->GetAtomOrthogonalCoordinate());
    return distance;
}

void PdbPreprocessor::ExtractCYSResidues(string pdb_file_path)
{
    PdbFile* pdb_file = new PdbFile(pdb_file_path);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();
    PdbFileSpace::PdbFile::PdbResidueVector cys_residues = GetAllCYSResidues(pdb_residues);
    PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = cys_residues.begin(); it != cys_residues.end(); it++)
    {
        PdbResidue* first_residue = (*it);
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = it+1 ; it1 != cys_residues.end(); it1++)
        {
            PdbResidue* second_residue = (*it1);
            double distance = GetDistanceofCYS(first_residue, second_residue, pdb_file, residue_atom_map);
            if (distance < dSulfurCutoff)
            {
                PdbPreprocessorDisulfideBond* disulfide_bond =
                        new PdbPreprocessorDisulfideBond(first_residue->GetResidueChainId(), second_residue->GetResidueChainId(),
                                                         first_residue->GetResidueSequenceNumber(), second_residue->GetResidueSequenceNumber(),
                                                         distance, true, first_residue->GetResidueInsertionCode(), second_residue->GetResidueInsertionCode(),
                                                         first_residue->GetResidueAlternateLocation(), second_residue->GetResidueAlternateLocation() );
                disulfide_bonds_.push_back(disulfide_bond);
            }
        }
    }
}

void PdbPreprocessor::UpdateCYSResidues(PdbFile *pdb_file, PdbPreprocessorDisulfideBondVector disulfide_bonds)
{
    for(PdbPreprocessorDisulfideBondVector::iterator it = disulfide_bonds.begin(); it != disulfide_bonds.end(); it++)
    {
        PdbPreprocessorDisulfideBond* disulfide_bond = (*it);
        if(disulfide_bond->GetIsBonded())
        {
            PdbAtom* pdb_atom_1 =
                    new PdbAtom(disulfide_bond->GetResidueChainId1(), "HG", "CYS",
                                disulfide_bond->GetResidueSequenceNumber1(), disulfide_bond->GetResidueInsertionCode1(), disulfide_bond->GetResidueAlternateLocation1());
            pdb_file->DeleteAtom(pdb_atom_1);
            PdbAtom* pdb_atom_2 =
                    new PdbAtom(disulfide_bond->GetResidueChainId2(), "HG", "CYS",
                                disulfide_bond->GetResidueSequenceNumber2(), disulfide_bond->GetResidueInsertionCode2(), disulfide_bond->GetResidueAlternateLocation2());
            pdb_file->DeleteAtom(pdb_atom_2);
            string target_key1;
            stringstream ss_1;
            ss_1 << "CYS" << "_" << disulfide_bond->GetResidueChainId1() << "_" << disulfide_bond->GetResidueSequenceNumber1() << disulfide_bond->GetResidueInsertionCode1()
                 << "_" << disulfide_bond->GetResidueAlternateLocation1();
            target_key1 = ss_1.str();
            string target_key2;
            stringstream ss_2;
            ss_2 << "CYS" << "_" << disulfide_bond->GetResidueChainId2() << "_" << disulfide_bond->GetResidueSequenceNumber2() << disulfide_bond->GetResidueInsertionCode2()
                 << "_" << disulfide_bond->GetResidueAlternateLocation2();
            target_key2 = ss_2.str();
            PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
            for(PdbFile::PdbResidueVector::iterator it1 = pdb_residues.begin(); it1 != pdb_residues.end(); it1++)
            {
                PdbResidue* pdb_residue = (*it1);

                string residue_name = pdb_residue->GetResidueName();
                if((residue_name).compare("CYS") == 0)
                {
                    stringstream ss;
                    string pdb_residue_key;
                    ss << residue_name << "_" << pdb_residue->GetResidueChainId() << "_" << pdb_residue->GetResidueSequenceNumber() << "_" << pdb_residue->GetResidueInsertionCode()
                          << pdb_residue->GetResidueAlternateLocation();
                    pdb_residue_key = ss.str();
                    if(pdb_residue_key.compare(target_key1) == 0 || pdb_residue_key.compare(target_key2) == 0)
                    {
                        //                        pdb_residue->SetResidueName("CYX");
                        pdb_file->UpdateResidueName(pdb_residue, "CYX");
                    }
                }
            }

        }
    }
}

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetAllHISResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues)
{
    PdbFileSpace::PdbFile::PdbResidueVector all_his_residues;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* pdb_residue = (*it);
        string pdb_residue_name = pdb_residue->GetResidueName();
        if((pdb_residue_name).compare("HIS") == 0)
        {
            all_his_residues.push_back(pdb_residue);
        }
    }
    return all_his_residues;
}

void PdbPreprocessor::ExtractHISResidues(string pdb_file_path)
{
    PdbFile* pdb_file = new PdbFile(pdb_file_path);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();
    PdbFileSpace::PdbFile::PdbResidueVector his_residues = GetAllHISResidues(pdb_residues);
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = his_residues.begin(); it != his_residues.end(); it++)
    {
        PdbResidue* his_residue = (*it);
        PdbPreprocessorHistidineMapping* histidine_mapping =
                new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HIE,
                                                    his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
        histidine_mappings_.push_back(histidine_mapping);
    }
}

void PdbPreprocessor::UpdateHISMapping(PdbFile *pdb_file, PdbPreprocessor::PdbPreprocessorHistidineMappingVector histidine_mappings)
{
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    for(PdbPreprocessorHistidineMappingVector::iterator it = histidine_mappings.begin(); it != histidine_mappings.end(); it++)
    {
        PdbPreprocessorHistidineMapping* histidine_mapping = (*it);
        string target_key;
        stringstream ss;
        ss << "HIS" << "_" << histidine_mapping->GetResidueChainId() << "_" << histidine_mapping->GetResidueSequenceNumber()
           << "_" << histidine_mapping->GetResidueInsertionCode() << "_" << histidine_mapping->GetResidueAlternateLocation();
        target_key = ss.str();
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = pdb_residues.begin(); it1 != pdb_residues.end(); it1++)
        {
            PdbResidue* pdb_residue = (*it1);

            string residue_name = pdb_residue->GetResidueName();
            if((residue_name).compare("HIS") == 0)
            {
                stringstream ss1;
                string pdb_residue_key;
                ss1 << residue_name << "_" << pdb_residue->GetResidueChainId() << "_" << pdb_residue->GetResidueSequenceNumber()
                    << "_" << pdb_residue->GetResidueInsertionCode() << "_" << pdb_residue->GetResidueAlternateLocation();
                pdb_residue_key = ss1.str();
                if(pdb_residue_key.compare(target_key) == 0)
                {
                    //                    pdb_residue->SetResidueName(histidine_mapping->GetStringFormatOfSelectedMapping());
                    pdb_file->UpdateResidueName(pdb_residue, histidine_mapping->GetStringFormatOfSelectedMapping());
                }
            }
        }
    }
}

vector<string> PdbPreprocessor::GetUnknownHeavyAtomNamesOfResidue(vector<string> pdb_atom_names_of_residue, vector<string> dataset_atom_names_of_residue)
{
    vector<string> unknown_heavy_atom_names_of_residue;
    for(vector<string>::iterator it = pdb_atom_names_of_residue.begin(); it != pdb_atom_names_of_residue.end(); it++)
    {
        string pdb_atom_name = (*it);
        bool found = false;
        if(!(pdb_atom_name.substr(0,1).compare("H") == 0 ||
             (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            for(vector<string>::iterator it1 = dataset_atom_names_of_residue.begin(); it1 != dataset_atom_names_of_residue.end(); it1++)
            {
                string dataset_atom_name = (*it1);
                {
                    if( pdb_atom_name.compare(dataset_atom_name) != 0)
                        continue;
                    else
                    {
                        found = true;
                        break;
                    }

                }
            }
            if(!found)
                unknown_heavy_atom_names_of_residue.push_back(pdb_atom_name);
        }
    }
    return unknown_heavy_atom_names_of_residue;
}
vector<string> PdbPreprocessor::GetAllAtomNamesOfResidueFromMultipleLibFiles(string residue_name, vector<string> lib_files)
{
    vector<string> all_atom_names_of_residue;
    bool found = false;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        if(!found)
        {
            LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
            vector<string> residue_names = lib_file->GetAllResidueNames();
            for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
            {
                string residue_name_from_lib = (*it);
                if((residue_name).compare(residue_name_from_lib) == 0)
                {
                    all_atom_names_of_residue = lib_file->GetAllAtomNamesOfResidue(residue_name);
                    found = true;
                    break;
                }
            }
        }
        else
            break;
    }
    return all_atom_names_of_residue;
}

vector<string> PdbPreprocessor::GetAllAtomNamesOfResidueFromMultiplePrepFiles(string residue_name, vector<string> prep_files)
{
    vector<string> all_atom_names_of_residue;
    bool found = false;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        if(!found)
        {
            PrepFileSpace::PrepFile* prep_files = new PrepFileSpace::PrepFile(*it);
            vector<string> residue_names = prep_files->GetAllResidueNames();
            for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
            {
                string residue_name_from_prep = (*it);
                if((residue_name).compare(residue_name_from_prep) == 0)
                {
                    all_atom_names_of_residue = prep_files->GetAllAtomNamesOfResidue(residue_name);
                    found = true;
                    break;
                }
            }
        }
        else
            break;
    }
    return all_atom_names_of_residue;
}

vector<string> PdbPreprocessor::GetAllAtomNamesOfResidueFromDatasetFiles(string residue_name, vector<string> lib_files, vector<string> prep_files)
{
    vector<string> all_atom_names_of_residue_from_lib = GetAllAtomNamesOfResidueFromMultipleLibFiles(residue_name, lib_files);
    vector<string> all_atom_names_of_residue_from_prep = GetAllAtomNamesOfResidueFromMultiplePrepFiles(residue_name, prep_files);
    vector<string> all_atom_names;
    if(all_atom_names_of_residue_from_prep.size() == 0)
    {
        for(vector<string>::iterator it1 = all_atom_names_of_residue_from_lib.begin(); it1 != all_atom_names_of_residue_from_lib.end(); it1++)
        {
            all_atom_names.push_back(*it1);
        }
    }
    else
    {
        for(vector<string>::iterator it1 = all_atom_names_of_residue_from_prep.begin(); it1 != all_atom_names_of_residue_from_prep.end(); it1++)
        {
            all_atom_names.push_back(*it1);
        }
    }
    return all_atom_names;
}

PdbFileSpace::PdbFile::PdbAtomVector PdbPreprocessor::GetUnknownHeavyAtomsOfResidue(PdbFile::PdbAtomVector pdb_atoms, vector<string> dataset_atom_names_of_residue)
{
    PdbFile::PdbAtomVector unknown_heavy_atoms_of_residue;
    for(PdbFile::PdbAtomVector::iterator it = pdb_atoms.begin(); it != pdb_atoms.end(); it++)
    {
        PdbAtom* pdb_atom = *it;
        string pdb_atom_name = pdb_atom->GetAtomName();
        bool found = false;
        if(!(pdb_atom_name.substr(0,1).compare("H") == 0 ||
             (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            for(vector<string>::iterator it1 = dataset_atom_names_of_residue.begin(); it1 != dataset_atom_names_of_residue.end(); it1++)
            {
                string dataset_atom_name = (*it1);

                if( pdb_atom_name.compare(dataset_atom_name) != 0)
                    continue;
                else
                {
                    found = true;
                    break;
                }
            }
            if(!found)
                unknown_heavy_atoms_of_residue.push_back(pdb_atom);
        }
    }
    return unknown_heavy_atoms_of_residue;
}

void PdbPreprocessor::ExtractUnknownHeavyAtoms(string pdb_file_path, vector<string> lib_files, vector<string> prep_files)
{
    PdbFile* pdb_file = new PdbFile(pdb_file_path);

    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);
cout << recognized_residues.size() << endl;
    PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* recognized_residue = *it;
        string residue_name = recognized_residue->GetResidueName();
        char chain_id = recognized_residue->GetResidueChainId();
        int sequence_number = recognized_residue->GetResidueSequenceNumber();
        char insertion_code = recognized_residue->GetResidueInsertionCode();
        char alternate_location = recognized_residue->GetResidueAlternateLocation();

        stringstream ss;
        ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
        string key = ss.str();
        PdbFile::PdbAtomVector atoms_of_residue = *(residue_atom_map[key]);
        vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
        PdbFile::PdbAtomVector unknown_heavy_atoms = GetUnknownHeavyAtomsOfResidue(atoms_of_residue, dataset_atom_names_of_residue);

        for(PdbFileSpace::PdbFile::PdbAtomVector::iterator it1 = unknown_heavy_atoms.begin(); it1 != unknown_heavy_atoms.end(); it1++)
        {
            PdbAtom* heavy_atom = (*it1);
            PdbPreprocessorUnrecognizedHeavyAtom* unknown_heavy_atom =
                    new PdbPreprocessorUnrecognizedHeavyAtom(heavy_atom->GetAtomChainId(), heavy_atom->GetAtomSerialNumber(),
                                                             heavy_atom->GetAtomName(), heavy_atom->GetAtomResidueName(),
                                                             heavy_atom->GetAtomResidueSequenceNumber(), heavy_atom->GetAtomInsertionCode(), heavy_atom->GetAtomAlternateLocation());
            unrecognized_heavy_atoms_.push_back(unknown_heavy_atom);
        }
    }
}

void PdbPreprocessor::RemoveUnknownHeavyAtoms(PdbFile *pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms)
{
    for(PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unknown_heavy_atoms.begin(); it != unknown_heavy_atoms.end(); it++)
    {
        PdbPreprocessorUnrecognizedHeavyAtom* unknown_heavy_atom = (*it);
        PdbAtom* pdb_atom =
                new PdbAtom(unknown_heavy_atom->GetResidueChainId(),unknown_heavy_atom->GetAtomName(),
                            unknown_heavy_atom->GetResidueName(), unknown_heavy_atom->GetResidueSequenceNumber(), unknown_heavy_atom->GetResidueInsertionCode(), unknown_heavy_atom->GetResidueAlternateLocation());
        pdb_file->DeleteAtom(pdb_atom);
    }
}

vector<string> PdbPreprocessor::GetRemovedHydrogenNamesOfResidue(vector<string> pdb_atom_names_of_residue, vector<string> dataset_atom_names_of_residue)
{
    vector<string> removed_hydrogen_names_of_residue;
    for(vector<string>::iterator it = pdb_atom_names_of_residue.begin(); it != pdb_atom_names_of_residue.end(); it++)
    {
        string pdb_atom_name = (*it);
        bool found = false;
        if((pdb_atom_name.substr(0,1).compare("H") == 0 ||
            (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            for(vector<string>::iterator it1 = dataset_atom_names_of_residue.begin(); it1 != dataset_atom_names_of_residue.end(); it1++)
            {
                string dataset_atom_name = (*it1);

                if( pdb_atom_name.compare(dataset_atom_name) != 0)
                    continue;
                else
                {
                    found = true;
                    break;
                }
            }

            if(!found)
                removed_hydrogen_names_of_residue.push_back(pdb_atom_name);
        }
    }
    return removed_hydrogen_names_of_residue;
}

PdbFileSpace::PdbFile::PdbAtomVector PdbPreprocessor::GetRemovedHydrogensOfResidue(PdbFile::PdbAtomVector pdb_atoms, vector<string> dataset_atom_names_of_residue)
{
    PdbFile::PdbAtomVector removed_hydrogens_of_residue;
    for(PdbFile::PdbAtomVector::iterator it = pdb_atoms.begin(); it != pdb_atoms.end(); it++)
    {
        PdbAtom* pdb_atom = *it;
        string pdb_atom_name = pdb_atom->GetAtomName();
        bool found = false;
        if((pdb_atom_name.substr(0,1).compare("H") == 0 ||
            (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            for(vector<string>::iterator it1 = dataset_atom_names_of_residue.begin(); it1 != dataset_atom_names_of_residue.end(); it1++)
            {
                string dataset_atom_name = (*it1);

                if( pdb_atom_name.compare(dataset_atom_name) != 0)
                    continue;
                else
                {
                    found = true;
                    break;
                }
            }
            if(!found)
                removed_hydrogens_of_residue.push_back(pdb_atom);
        }
    }
    return removed_hydrogens_of_residue;
}

void PdbPreprocessor::ExtractRemovedHydrogens(string pdb_file_path, vector<string> lib_files, vector<string> prep_files)
{
    PdbFile* pdb_file = new PdbFile(pdb_file_path);

    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

    PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* recognized_residue = *it;
        string residue_name = recognized_residue->GetResidueName();
        char chain_id = recognized_residue->GetResidueChainId();
        int sequence_number = recognized_residue->GetResidueSequenceNumber();
        char insertion_code = recognized_residue->GetResidueInsertionCode();
        char alternate_location = recognized_residue->GetResidueAlternateLocation();
        stringstream ss;
        ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
        string key = ss.str();
        PdbFile::PdbAtomVector atoms_of_residue = *(residue_atom_map[key]);
        vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
        PdbFile::PdbAtomVector removed_hydrogens = GetRemovedHydrogensOfResidue(atoms_of_residue, dataset_atom_names_of_residue);
        for(PdbFileSpace::PdbFile::PdbAtomVector::iterator it1 = removed_hydrogens.begin(); it1 != removed_hydrogens.end(); it1++)
        {
            PdbAtom* removed_hydrogen = (*it1);
            PdbPreprocessorReplacedHydrogen* removed_hydrogen_atom =
                    new PdbPreprocessorReplacedHydrogen(removed_hydrogen->GetAtomChainId(), removed_hydrogen->GetAtomSerialNumber(), removed_hydrogen->GetAtomName(),
                                                        removed_hydrogen->GetAtomResidueName(), removed_hydrogen->GetAtomResidueSequenceNumber(),
                                                        removed_hydrogen->GetAtomInsertionCode(), removed_hydrogen->GetAtomAlternateLocation());
            replaced_hydrogens_.push_back(removed_hydrogen_atom);
        }
    }
}

void PdbPreprocessor::RemoveRemovedHydrogens(PdbFile *pdb_file, PdbPreprocessorReplacedHydrogenVector replaced_hydrogens)
{
    for(PdbPreprocessorReplacedHydrogenVector::iterator it = replaced_hydrogens.begin(); it != replaced_hydrogens.end(); it++)
    {
        PdbPreprocessorReplacedHydrogen* replaced_hydrogen = (*it);
        PdbAtom* pdb_atom =
                new PdbAtom(replaced_hydrogen->GetResidueChainId(),replaced_hydrogen->GetAtomName(),
                            replaced_hydrogen->GetResidueName(), replaced_hydrogen->GetResidueSequenceNumber(), replaced_hydrogen->GetResidueInsertionCode(), replaced_hydrogen->GetResidueAlternateLocation());
        pdb_file->DeleteAtom(pdb_atom);
    }
}

void PdbPreprocessor::ExtractAminoAcidChains(string pdb_file_path)
{
    PdbFile* pdb_file = new PdbFile(pdb_file_path);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();

    PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
    PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* residue = *it;
        char chain_id = residue->GetResidueChainId();
        int sequence_number = residue->GetResidueSequenceNumber();
        char insertion_code = residue->GetResidueInsertionCode();

        chain_map_sequence_number[chain_id].push_back(sequence_number);
        chain_map_insertion_code[chain_id].push_back(insertion_code);
    }

    for(PdbPreprocessorChainIdSequenceNumbersMap::iterator it = chain_map_sequence_number.begin(); it != chain_map_sequence_number.end(); it++)
    {
        char chain_id = (*it).first;
        vector<int> sequence_numbers = (*it).second;
        vector<char> insertion_codes = chain_map_insertion_code[chain_id];

        vector<int>::iterator starting_sequence_number_iterator = min_element(sequence_numbers.begin(), sequence_numbers.end());
        vector<int>::iterator ending_sequence_number_iterator = max_element(sequence_numbers.begin(), sequence_numbers.end());
        int starting_index = distance(sequence_numbers.begin(), starting_sequence_number_iterator);
        int ending_index = distance(sequence_numbers.begin(), ending_sequence_number_iterator);

        PdbPreprocessorChainTermination* chain = new PdbPreprocessorChainTermination(chain_id, *starting_sequence_number_iterator, *ending_sequence_number_iterator,
                                                                                     insertion_codes.at(starting_index), insertion_codes.at(ending_index) );

        chain_terminations_.push_back(chain);
    }
}

void PdbPreprocessor::UpdateAminoAcidChains(PdbFile *pdb_file, vector<string> lib_files, PdbPreprocessorChainTerminationVector chain_terminations)
{
    for(PdbPreprocessor::PdbPreprocessorChainTerminationVector::iterator it1 = chain_terminations.begin(); it1 != chain_terminations.end(); it1++)
    {
        PdbPreprocessorChainTermination* chain = (*it1);
        // Zwitterionic in n terminal
        if(chain->GetStringFormatOfSelectedNTermination().find("+") != string::npos || chain->GetStringFormatOfSelectedNTermination().find("-") != string::npos)
        {
            // Zwitterionic in c terminal
            if(chain->GetStringFormatOfSelectedCTermination().find("+") != string::npos || chain->GetStringFormatOfSelectedCTermination().find("-") != string::npos)
            {
                // End of chain
                // Do nothing
            }
            else
            {
                // Add c terminal at the end of the chain
                PossibleCChainTermination c_termination = chain->GetSelectedCTermination();
                string string_c_termination = chain->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, chain->GetResidueChainId(),
                                                                      chain->GetEndingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_file->InsertResidueAfter(pdb_atom_card);
                }
            }
        }
        else
        {
            // Zwitterionic in c terminal
            if(chain->GetStringFormatOfSelectedCTermination().find("+") != string::npos || chain->GetStringFormatOfSelectedCTermination().find("-") != string::npos)
            {
                // Add n terminal residue at the beginning of the chain
                PossibleNChainTermination n_termination = chain->GetSelectedNTermination();
                string string_n_termination = chain->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, chain->GetResidueChainId(),
                                                                      chain->GetStartingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_file->InsertResidueBefore(pdb_atom_card);
                }
            }
            else
            {
                // Add c terminal residue at the end of the chain
                PossibleCChainTermination c_termination = chain->GetSelectedCTermination();
                string string_c_termination = chain->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue_from_c_termination = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, lib_files);
                if(lib_file_residue_from_c_termination != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms_from_c_termination = lib_file_residue_from_c_termination->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card_for_c_termination = new PdbAtomCard();
                    pdb_atom_card_for_c_termination->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map_for_c_termination;
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it3 = lib_atoms_from_c_termination.begin(); it3 != lib_atoms_from_c_termination.end(); it3++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it3).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, chain->GetResidueChainId(),
                                                                      chain->GetEndingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map_for_c_termination[serial_number] = pdb_atom;
                        serial_number++;
                    }
                    pdb_atom_card_for_c_termination->SetAtoms(atom_map_for_c_termination);
                    pdb_file->InsertResidueAfter(pdb_atom_card_for_c_termination);
                }
                // Add n terminal residue at the beginning of the chain
                PossibleNChainTermination n_termination = chain->GetSelectedNTermination();
                string string_n_termination = chain->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, chain->GetResidueChainId(),
                                                                      chain->GetStartingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_file->InsertResidueBefore(pdb_atom_card);
                }
            }
        }
    }
}

void PdbPreprocessor::ExtractGapsInAminoAcidChains(string pdb_file_path)
{
    PdbFile* pdb_file = new PdbFile(pdb_file_path);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();

    PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
    PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* residue = *it;
        char chain_id = residue->GetResidueChainId();
        int sequence_number = residue->GetResidueSequenceNumber();
        char insertion_code = residue->GetResidueInsertionCode();

        chain_map_sequence_number[chain_id].push_back(sequence_number);
        chain_map_insertion_code[chain_id].push_back(insertion_code);
    }

    for(PdbPreprocessorChainIdSequenceNumbersMap::iterator it = chain_map_sequence_number.begin(); it != chain_map_sequence_number.end(); it++)
    {
        char chain_id = (*it).first;
        vector<int> sequence_numbers = (*it).second;
        vector<char> insertion_codes = chain_map_insertion_code[chain_id];

        int starting_sequence_number = *min_element(sequence_numbers.begin(), sequence_numbers.end());
        int ending_sequence_number = *max_element(sequence_numbers.begin(), sequence_numbers.end());

        for(unsigned int i = 0; i < sequence_numbers.size() - 1; i++)
        {
            if(sequence_numbers.at(i) == sequence_numbers.at(i+1) || (sequence_numbers.at(i) + 1) == sequence_numbers.at(i+1))
            {

            }
            else
            {
                PdbPreprocessorMissingResidue* missing_residues = new PdbPreprocessorMissingResidue(chain_id, starting_sequence_number,
                                                                                                    ending_sequence_number, sequence_numbers.at(i),
                                                                                                    sequence_numbers.at(i+1), insertion_codes.at(i), insertion_codes.at(i+1));
                missing_residues_.push_back(missing_residues);
            }
        }
    }
}

void PdbPreprocessor::UpdateGapsInAminoAcidChains(PdbFile* pdb_file, vector<string> lib_files, PdbPreprocessorMissingResidueVector gaps)
{
    for(PdbPreprocessor::PdbPreprocessorMissingResidueVector::iterator it1 = gaps.begin(); it1 != gaps.end(); it1++)
    {
        PdbPreprocessorMissingResidue* gap = (*it1);
        pdb_file->SplitAtomCardOfModelCard(gap->GetResidueChainId(), gap->GetResidueAfterGap());
    }
    for(PdbPreprocessor::PdbPreprocessorMissingResidueVector::iterator it1 = gaps.begin(); it1 != gaps.end(); it1++)
    {
        PdbPreprocessorMissingResidue* gap = (*it1);
        // Zwitterionic in n terminal
        if(gap->GetStringFormatOfSelectedNTermination().find("+") != string::npos || gap->GetStringFormatOfSelectedNTermination().find("-") != string::npos)
        {
            // Zwitterionic in c terminal
            if(gap->GetStringFormatOfSelectedCTermination().find("+") != string::npos || gap->GetStringFormatOfSelectedCTermination().find("-") != string::npos)
            {

            }
            else
            {
                // Add c terminal at the end of the chain
                PossibleCChainTermination c_termination = gap->GetSelectedCTermination();
                string string_c_termination = gap->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueBeforeGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_file->InsertResidueAfter(pdb_atom_card);
                }
            }
        }
        else
        {
            // Zwitterionic in c terminal
            if(gap->GetStringFormatOfSelectedCTermination().find("+") != string::npos || gap->GetStringFormatOfSelectedCTermination().find("-") != string::npos)
            {
                // Add n terminal residue at the beginning of the chain
                PossibleNChainTermination n_termination = gap->GetSelectedNTermination();
                string string_n_termination = gap->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueAfterGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_file->InsertResidueBefore(pdb_atom_card);
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
                // Add c terminal residue at the end of the chain
                PossibleCChainTermination c_termination = gap->GetSelectedCTermination();
                string string_c_termination = gap->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue_from_c_termination = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, lib_files);
                if(lib_file_residue_from_c_termination != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms_from_c_termination = lib_file_residue_from_c_termination->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card_for_c_termination = new PdbAtomCard();
                    pdb_atom_card_for_c_termination->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map_for_c_termination;
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it3 = lib_atoms_from_c_termination.begin(); it3 != lib_atoms_from_c_termination.end(); it3++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it3).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueBeforeGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map_for_c_termination[serial_number] = pdb_atom;
                        serial_number++;
                    }
                    pdb_atom_card_for_c_termination->SetAtoms(atom_map_for_c_termination);
                    pdb_file->InsertResidueAfter(pdb_atom_card_for_c_termination);
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
                // Add n terminal residue at the beginning of the chain
                PossibleNChainTermination n_termination = gap->GetSelectedNTermination();
                string string_n_termination = gap->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueAfterGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_file->InsertResidueBefore(pdb_atom_card);
                }
            }
        }
    }

}

LibraryFileSpace::LibraryFileResidue* PdbPreprocessor::GetLibraryResidueByNameFromMultipleLibraryFiles(string residue_name, vector<string> lib_files)
{
    LibraryFileSpace::LibraryFileResidue* library_residue = new LibraryFileSpace::LibraryFileResidue();
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        library_residue = lib_file->GetLibraryResidueByResidueName(residue_name);
        if(library_residue != NULL)
            break;
    }
    return library_residue;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessor::Print(ostream &out)
{
}

