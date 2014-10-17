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
#include "../../../includes/FileSet/PdbFileSpace/pdbresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/common.hpp"
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
PdbPreprocessor::PdbPreprocessorAlternateResidueMap PdbPreprocessor::GetAlternateResidueMap(){
    return alternate_residue_map_;
}
PdbPreprocessor::PdbPreprocessorToBeDeletedAtomVector PdbPreprocessor::GetToBeDeletedAtoms(){
    return to_be_deleted_atoms_;
}
PdbPreprocessor::PdbPreprocessorToBeDeletedResidueVector PdbPreprocessor::GetToBeDeletedResidues(){
    return to_be_deleted_residues_;
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
void PdbPreprocessor::SetToBeDeletedAtoms(PdbPreprocessorToBeDeletedAtomVector to_be_deleted_atoms)
{
    to_be_deleted_atoms_.clear();
    for(PdbPreprocessorToBeDeletedAtomVector::iterator it = to_be_deleted_atoms.begin(); it != to_be_deleted_atoms.end(); it++)
    {
        to_be_deleted_atoms_.push_back(*it);
    }
}
void PdbPreprocessor::SetToBeDeletedResidues(PdbPreprocessorToBeDeletedResidueVector to_be_deleted_residues)
{
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorToBeDeletedResidueVector::iterator it = to_be_deleted_residues.begin(); it != to_be_deleted_residues.end(); it++)
    {
        to_be_deleted_residues_.push_back(*it);
    }
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
vector<string> PdbPreprocessor::GetUnrecognizedResidueNames(vector<string> pdb_residue_names, vector<string> dataset_residue_names)
{
    vector<string> unrecognized_residue_names;
    bool is_recognized = false;
    int counter = 0;
    for(vector<string>::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        string pdb_residue_name = *it;
        if(pdb_residue_name.compare("HIS") != 0 )
        {
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
        else
            counter++;
    }
    cout << counter << " HIS residue(s) found" << endl;
    return unrecognized_residue_names;
}

ResidueNameMap PdbPreprocessor::GetUnrecognizedResidueNamesMap(vector<string> pdb_residue_names, ResidueNameMap dataset_residue_names)
{
    ResidueNameMap unrecognized_residue_names = ResidueNameMap();
    int counter = 0;
    for(vector<string>::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        string pdb_residue_name = *it;
        if(pdb_residue_name.compare("HIS") != 0 )
        {
            if(dataset_residue_names.find(pdb_residue_name) == dataset_residue_names.end())
            {
                unrecognized_residue_names[pdb_residue_name] = pdb_residue_name;
            }
        }
        else
            counter++;
    }
    cout << counter << " HIS residue(s) found" << endl;
    return unrecognized_residue_names;
}

vector<string> PdbPreprocessor::GetRecognizedResidueNames(vector<string> pdb_residue_names, vector<string> dataset_residue_names)
{
    vector<string> recognized_residue_names;
    int counter = 0;
    for(vector<string>::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        string pdb_residue_name = *it;
        if(pdb_residue_name.compare("HIS") != 0)
        {
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
        else
        {
            recognized_residue_names.push_back(pdb_residue_name);
            counter++;
        }
    }
    cout << counter << " HIS residue(s) found" << endl;
    return recognized_residue_names;
}

ResidueNameMap PdbPreprocessor::GetRecognizedResidueNamesMap(vector<string> pdb_residue_names, ResidueNameMap dataset_residue_names)
{
    ResidueNameMap recognized_residue_names = ResidueNameMap();
    int counter = 0;
    for(vector<string>::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        string pdb_residue_name = *it;
        if(pdb_residue_name.compare("HIS") != 0)
        {
            if(dataset_residue_names.find(pdb_residue_name) != dataset_residue_names.end())
            {
                recognized_residue_names[pdb_residue_name] = pdb_residue_name;
            }
        }
        else
        {
            recognized_residue_names[pdb_residue_name] = pdb_residue_name;
            counter++;
        }
    }
    cout << counter << " HIS residue(s) found" << endl;
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

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetUnrecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, ResidueNameMap unrecognized_residue_names)
{
    PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* pdb_residue = *it;
        string pdb_residue_name = pdb_residue->GetResidueName();
        if(unrecognized_residue_names.find(pdb_residue_name) != unrecognized_residue_names.end())
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
        {
            recognized_residues.push_back(pdb_residue);
        }
    }
    return recognized_residues;
}

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetRecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, ResidueNameMap recognized_residue_names)
{
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* pdb_residue = *it;
        string pdb_residue_name = pdb_residue->GetResidueName();
        if(recognized_residue_names.find(pdb_residue_name) != recognized_residue_names.end())
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

ResidueNameMap PdbPreprocessor::GetAllResidueNamesFromMultipleLibFilesMap(vector<string> lib_files)
{
    ResidueNameMap all_residue_names;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        vector<string> residue_names = lib_file->GetAllResidueNames();
        for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
        {
            string residue_name = (*it1);
            all_residue_names[residue_name] = (residue_name);
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

ResidueNameMap PdbPreprocessor::GetAllResidueNamesFromMultiplePrepFilesMap(vector<string> prep_files)
{
    ResidueNameMap all_residue_names;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        vector<string> residue_names = prep_file->GetAllResidueNames();
        for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
        {
            string residue_name = (*it1);
            all_residue_names[residue_name] = residue_name;
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

ResidueNameMap PdbPreprocessor::GetAllResidueNamesFromDatasetFilesMap(vector<string> lib_files, vector<string> prep_files)
{
    ResidueNameMap all_lib_residue_names = GetAllResidueNamesFromMultipleLibFilesMap(lib_files);
    ResidueNameMap all_prep_residue_names = GetAllResidueNamesFromMultiplePrepFilesMap(prep_files);
    ResidueNameMap all_residue_names = ResidueNameMap();
    for(ResidueNameMap::iterator it1 = all_lib_residue_names.begin(); it1 != all_lib_residue_names.end(); it1++)
    {
        string key = (*it1).first;
        string val = (*it1).second;
        all_residue_names[key] = val;
    }
    for(ResidueNameMap::iterator it1 = all_prep_residue_names.begin(); it1 != all_prep_residue_names.end(); it1++)
    {
        string key = (*it1).first;
        string val = (*it1).second;
        all_residue_names[key] = val;
    }
    return all_residue_names;
}

void PdbPreprocessor::ExtractUnrecognizedResidues(string pdb_file_path, vector<string> lib_files, vector<string> prep_files)
{
    try
    {
        // Slow version
//        vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
        // Advanced version
        ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
        // Slow version
//        vector<string> unrecognized_residue_names = GetUnrecognizedResidueNames(pdb_residue_names, dataset_residue_names);
        // Advanced version
        ResidueNameMap unrecognized_residue_names = GetUnrecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);

        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues = GetUnrecognizedResidues(pdb_residues, unrecognized_residue_names);
        PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
        {
            PdbFileSpace::PdbResidue* residue = *it;
            char chain_id = residue->GetResidueChainId();
            int sequence_number = residue->GetResidueSequenceNumber();
            chain_map_sequence_number[chain_id].push_back(sequence_number);
        }
        PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_start_end_sequence_number;
        for(PdbPreprocessorChainIdSequenceNumbersMap::iterator it = chain_map_sequence_number.begin(); it != chain_map_sequence_number.end(); it++)
        {
            char chain_id = (*it).first;
            vector<int> sequence_numbers = (*it).second;
            chain_map_start_end_sequence_number[chain_id] = vector<int>();
            int starting_sequence_number = *min_element(sequence_numbers.begin(), sequence_numbers.end());
            int ending_sequence_number = *max_element(sequence_numbers.begin(), sequence_numbers.end());
            chain_map_start_end_sequence_number[chain_id].push_back(starting_sequence_number);
            chain_map_start_end_sequence_number[chain_id].push_back(ending_sequence_number);
        }
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
        {
            PdbResidue* pdb_residue = (*it);
            vector<int> residue_start_end_squence_number = chain_map_start_end_sequence_number[pdb_residue->GetResidueChainId()];
            PdbPreprocessorUnrecognizedResidue* unrecognized_residue =
                    new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(), pdb_residue->GetResidueSequenceNumber(),
                                                           pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation(), false);
            if(pdb_residue->GetResidueSequenceNumber() > residue_start_end_squence_number[0] && pdb_residue->GetResidueSequenceNumber() < residue_start_end_squence_number[1])
                unrecognized_residue->SetMiddleOfChain(true);
            unrecognized_residues_.push_back(unrecognized_residue);
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractUnrecognizedResidues(PdbFile* pdb_file, vector<string> lib_files, vector<string> prep_files)
{
    // Slow version
//    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    // Advanced version
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    // Slow version
//    vector<string> unrecognized_residue_names = GetUnrecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    // Advanced version
    ResidueNameMap unrecognized_residue_names = GetUnrecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);

    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues = GetUnrecognizedResidues(pdb_residues, unrecognized_residue_names);
    PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* residue = *it;
        char chain_id = residue->GetResidueChainId();
        int sequence_number = residue->GetResidueSequenceNumber();
        chain_map_sequence_number[chain_id].push_back(sequence_number);
    }
    PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_start_end_sequence_number;
    for(PdbPreprocessorChainIdSequenceNumbersMap::iterator it = chain_map_sequence_number.begin(); it != chain_map_sequence_number.end(); it++)
    {
        char chain_id = (*it).first;
        vector<int> sequence_numbers = (*it).second;
        chain_map_start_end_sequence_number[chain_id] = vector<int>();
        int starting_sequence_number = *min_element(sequence_numbers.begin(), sequence_numbers.end());
        int ending_sequence_number = *max_element(sequence_numbers.begin(), sequence_numbers.end());
        chain_map_start_end_sequence_number[chain_id].push_back(starting_sequence_number);
        chain_map_start_end_sequence_number[chain_id].push_back(ending_sequence_number);
    }
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
    {
        PdbResidue* pdb_residue = (*it);
        vector<int> residue_start_end_squence_number = chain_map_start_end_sequence_number[pdb_residue->GetResidueChainId()];
        PdbPreprocessorUnrecognizedResidue* unrecognized_residue =
                new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(), pdb_residue->GetResidueSequenceNumber(),
                                                       pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation(), false);
        if(pdb_residue->GetResidueSequenceNumber() > residue_start_end_squence_number[0] && pdb_residue->GetResidueSequenceNumber() < residue_start_end_squence_number[1])
            unrecognized_residue->SetMiddleOfChain(true);
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
    try
    {
        // Slow version
//        vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
        // Advanced version
        ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
        // Slow version
//        vector<string> unrecognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
        // Advanced version
        ResidueNameMap unrecognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);

        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, unrecognized_residue_names);

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
        {
            PdbResidue* pdb_residue = (*it);
            PdbPreprocessorUnrecognizedResidue* recognized_residue =
                    new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(), pdb_residue->GetResidueSequenceNumber(),
                                                           pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation(), false);
            recognized_residues_.push_back(recognized_residue);
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractRecognizedResidues(PdbFile* pdb_file, vector<string> lib_files, vector<string> prep_files)
{
    // Slow version
//    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    // Advanced version
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    // Slow version
//    vector<string> unrecognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    // Advanced version
    ResidueNameMap unrecognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, unrecognized_residue_names);

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        PdbResidue* pdb_residue = (*it);
        PdbPreprocessorUnrecognizedResidue* recognized_residue =
                new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(), pdb_residue->GetResidueSequenceNumber(),
                                                       pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation(), false);
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
    try
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
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractCYSResidues(PdbFile* pdb_file)
{
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
            ss_1 << "CYS" << "_" << disulfide_bond->GetResidueChainId1() << "_" << disulfide_bond->GetResidueSequenceNumber1() << "_" << disulfide_bond->GetResidueInsertionCode1()
                 << "_" << disulfide_bond->GetResidueAlternateLocation1();
            target_key1 = ss_1.str();
            string target_key2;
            stringstream ss_2;
            ss_2 << "CYS" << "_" << disulfide_bond->GetResidueChainId2() << "_" << disulfide_bond->GetResidueSequenceNumber2() << "_" << disulfide_bond->GetResidueInsertionCode2()
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
                       << "_" << pdb_residue->GetResidueAlternateLocation();
                    pdb_residue_key = ss.str();
                    if(pdb_residue_key.compare(target_key1) == 0 || pdb_residue_key.compare(target_key2) == 0)
                    {
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
    try
    {
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();
        PdbFileSpace::PdbFile::PdbResidueVector his_residues = GetAllHISResidues(pdb_residues);
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = his_residues.begin(); it != his_residues.end(); it++)
        {
            PdbResidue* his_residue = (*it);

            if(pdb_file->GetAtomOfResidueByName(his_residue, "HE2") != NULL && pdb_file->GetAtomOfResidueByName(his_residue, "HD1") == NULL)
            {
                PdbPreprocessorHistidineMapping* histidine_mapping =
                        new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HIE,
                                                            his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
                histidine_mappings_.push_back(histidine_mapping);
            }
            else if(pdb_file->GetAtomOfResidueByName(his_residue, "HE2") == NULL && pdb_file->GetAtomOfResidueByName(his_residue, "HD1") != NULL)
            {
                PdbPreprocessorHistidineMapping* histidine_mapping =
                        new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HID,
                                                            his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
                histidine_mappings_.push_back(histidine_mapping);
            }
            else if(pdb_file->GetAtomOfResidueByName(his_residue, "HE2") != NULL && pdb_file->GetAtomOfResidueByName(his_residue, "HD1") != NULL)
            {
                PdbPreprocessorHistidineMapping* histidine_mapping =
                        new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HIP,
                                                            his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
                histidine_mappings_.push_back(histidine_mapping);
            }
            else if(pdb_file->GetAtomOfResidueByName(his_residue, "HE2") == NULL && pdb_file->GetAtomOfResidueByName(his_residue, "HD1") == NULL)
            {
                PdbPreprocessorHistidineMapping* histidine_mapping =
                        new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HIE,
                                                            his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
                histidine_mappings_.push_back(histidine_mapping);
            }
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractHISResidues(PdbFile* pdb_file)
{
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();
    PdbFileSpace::PdbFile::PdbResidueVector his_residues = GetAllHISResidues(pdb_residues);
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = his_residues.begin(); it != his_residues.end(); it++)
    {
        PdbResidue* his_residue = (*it);

        // HIE residue
        if(pdb_file->GetAtomOfResidueByName(his_residue, "HE2") != NULL && pdb_file->GetAtomOfResidueByName(his_residue, "HD1") == NULL)
        {
            PdbPreprocessorHistidineMapping* histidine_mapping =
                    new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HIE,
                                                        his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
            histidine_mappings_.push_back(histidine_mapping);
        }
        // HID residue
        else if(pdb_file->GetAtomOfResidueByName(his_residue, "HE2") == NULL && pdb_file->GetAtomOfResidueByName(his_residue, "HD1") != NULL)
        {
            PdbPreprocessorHistidineMapping* histidine_mapping =
                    new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HID,
                                                        his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
            histidine_mappings_.push_back(histidine_mapping);
        }
        // HIP residue
        else if(pdb_file->GetAtomOfResidueByName(his_residue, "HE2") != NULL && pdb_file->GetAtomOfResidueByName(his_residue, "HD1") != NULL)
        {
            PdbPreprocessorHistidineMapping* histidine_mapping =
                    new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HIP,
                                                        his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
            histidine_mappings_.push_back(histidine_mapping);
        }
        else if(pdb_file->GetAtomOfResidueByName(his_residue, "HE2") == NULL && pdb_file->GetAtomOfResidueByName(his_residue, "HD1") == NULL)
        {
            PdbPreprocessorHistidineMapping* histidine_mapping =
                    new PdbPreprocessorHistidineMapping(his_residue->GetResidueChainId(), his_residue->GetResidueSequenceNumber(), HIE,
                                                        his_residue->GetResidueInsertionCode(), his_residue->GetResidueAlternateLocation());
            histidine_mappings_.push_back(histidine_mapping);
        }
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
                    // HIE residue
                    PdbAtom* HE2 = pdb_file->GetAtomOfResidueByName(pdb_residue, "HE2");
                    PdbAtom* HD1 = pdb_file->GetAtomOfResidueByName(pdb_residue, "HD1");
                    if(HE2 != NULL && HD1 == NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == HID)
                        {
                            // Delete HE2
                            pdb_file->DeleteAtom(HE2);
                        }
                    }
                    // HID residue
                    if(HE2 == NULL && HD1 != NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == HIE)
                        {
                            // Delete HD1
                            pdb_file->DeleteAtom(HD1);
                        }

                    }
                    // HIP residue
                    if(HE2 != NULL && HD1 != NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == HIE)
                        {
                            // Delete HD1
                            pdb_file->DeleteAtom(HD1);
                        }
                        if(histidine_mapping->GetSelectedMapping() == HID)
                        {
                            // Delete HE2
                            pdb_file->DeleteAtom(HE2);
                        }
                    }
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

AtomNameMap PdbPreprocessor::GetUnknownHeavyAtomNamesOfResidue(vector<string> pdb_atom_names_of_residue, AtomNameMap dataset_atom_names_of_residue)
{
    AtomNameMap unknown_heavy_atom_names_of_residue = AtomNameMap();
    for(vector<string>::iterator it = pdb_atom_names_of_residue.begin(); it != pdb_atom_names_of_residue.end(); it++)
    {
        string pdb_atom_name = (*it);
        if(!(pdb_atom_name.substr(0,1).compare("H") == 0 ||
             (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(dataset_atom_names_of_residue.find(pdb_atom_name) == dataset_atom_names_of_residue.end())
                unknown_heavy_atom_names_of_residue[pdb_atom_name] = pdb_atom_name;
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
                string residue_name_from_lib = (*it1);
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
//
AtomNameMap PdbPreprocessor::GetAllAtomNamesOfResidueFromMultipleLibFilesMap(string residue_name, vector<string> lib_files)
{
    AtomNameMap all_atom_names_of_residue = AtomNameMap();
    bool found = false;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        if(!found)
        {
            LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
            ResidueNameMap residue_names = lib_file->GetAllResidueNamesMap();
            if(residue_names.find(residue_name) != residue_names.end())
            {
                all_atom_names_of_residue = lib_file->GetAllAtomNamesOfResidueMap(residue_name);
                found = true;
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
            PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
            vector<string> residue_names = prep_file->GetAllResidueNames();
            for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
            {
                string residue_name_from_prep = (*it1);
                if((residue_name).compare(residue_name_from_prep) == 0)
                {
                    all_atom_names_of_residue = prep_file->GetAllAtomNamesOfResidue(residue_name);
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
//
AtomNameMap PdbPreprocessor::GetAllAtomNamesOfResidueFromMultiplePrepFilesMap(string residue_name, vector<string> prep_files)
{
    AtomNameMap all_atom_names_of_residue = AtomNameMap();
    bool found = false;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        if(!found)
        {
            PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
            ResidueNameMap residue_names = prep_file->GetAllResidueNamesMap();
            if(residue_names.find(residue_name) != residue_names.end())
            {
                all_atom_names_of_residue = prep_file->GetAllAtomNamesOfResidueMap(residue_name);
                found = true;
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
    for(vector<string>::iterator it1 = all_atom_names_of_residue_from_lib.begin(); it1 != all_atom_names_of_residue_from_lib.end(); it1++)
    {
        all_atom_names.push_back(*it1);
    }
    for(vector<string>::iterator it1 = all_atom_names_of_residue_from_prep.begin(); it1 != all_atom_names_of_residue_from_prep.end(); it1++)
    {
        all_atom_names.push_back(*it1);
    }
    return all_atom_names;
}

AtomNameMap PdbPreprocessor::GetAllAtomNamesOfResidueFromDatasetFilesMap(string residue_name, vector<string> lib_files, vector<string> prep_files)
{
    AtomNameMap all_atom_names_of_residue_from_lib = GetAllAtomNamesOfResidueFromMultipleLibFilesMap(residue_name, lib_files);
    AtomNameMap all_atom_names_of_residue_from_prep = GetAllAtomNamesOfResidueFromMultiplePrepFilesMap(residue_name, prep_files);
    AtomNameMap all_atom_names = AtomNameMap();
    for(AtomNameMap::iterator it1 = all_atom_names_of_residue_from_lib.begin(); it1 != all_atom_names_of_residue_from_lib.end(); it1++)
    {
        string key = (*it1).first;
        string val = (*it1).second;
        all_atom_names[key] = val;
    }
    for(AtomNameMap::iterator it1 = all_atom_names_of_residue_from_prep.begin(); it1 != all_atom_names_of_residue_from_prep.end(); it1++)
    {
        string key = (*it1).first;
        string val = (*it1).second;
        all_atom_names[key] = val;
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
        if(!(pdb_atom_name.substr(0,1).compare("H") == 0 ||
             (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(!(find(dataset_atom_names_of_residue.begin(), dataset_atom_names_of_residue.end(), pdb_atom_name) != dataset_atom_names_of_residue.end()))
            {
                unknown_heavy_atoms_of_residue.push_back(pdb_atom);
            }
        }
    }
    return unknown_heavy_atoms_of_residue;
}

PdbFileSpace::PdbFile::PdbAtomVector PdbPreprocessor::GetUnknownHeavyAtomsOfResidue(PdbFile::PdbAtomVector pdb_atoms, AtomNameMap dataset_atom_names_of_residue)
{
    PdbFile::PdbAtomVector unknown_heavy_atoms_of_residue;
    for(PdbFile::PdbAtomVector::iterator it = pdb_atoms.begin(); it != pdb_atoms.end(); it++)
    {
        PdbAtom* pdb_atom = *it;
        string pdb_atom_name = pdb_atom->GetAtomName();
        if(!(pdb_atom_name.substr(0,1).compare("H") == 0 ||
             (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(dataset_atom_names_of_residue.find(pdb_atom_name) == dataset_atom_names_of_residue.end())
            {
                unknown_heavy_atoms_of_residue.push_back(pdb_atom);
            }
        }
    }
    return unknown_heavy_atoms_of_residue;
}

void PdbPreprocessor::ExtractUnknownHeavyAtoms(string pdb_file_path, vector<string> lib_files, vector<string> prep_files)
{
    try
    {
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        // Slow version
//         vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
        // Advanced version
        ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
        vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
        // Slow version
//         vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
        // Advanced version
        ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);


        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

        if(recognized_residues.size() > PdbResidueThreshold)
        {
            cout << "Number of recognized residues in this pdb file is more than a threshold." << endl;
            cout << "If you are using prep file, this may take a while, please wait ..." << endl;
        }

        PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
        {
            PdbFileSpace::PdbResidue* recognized_residue = *it;
            string residue_name = recognized_residue->GetResidueName();
            if(residue_name.compare("HIS") != 0)
            {
                char chain_id = recognized_residue->GetResidueChainId();
                int sequence_number = recognized_residue->GetResidueSequenceNumber();
                char insertion_code = recognized_residue->GetResidueInsertionCode();
                char alternate_location = recognized_residue->GetResidueAlternateLocation();

                stringstream ss;
                ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                string key = ss.str();
                PdbFile::PdbAtomVector atoms_of_residue = *(residue_atom_map[key]);
                // Slow version
//                vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
                // Advanced version
                AtomNameMap dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFilesMap(residue_name, lib_files, prep_files);


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
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractUnknownHeavyAtoms(PdbFile* pdb_file, vector<string> lib_files, vector<string> prep_files)
{
    // Slow version
//    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    // Advanced version
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    // Slow version
//    vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    // Advanced version
    ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);


    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

    if(recognized_residues.size() > PdbResidueThreshold)
    {
        cout << "Number of recognized residues in this pdb file is more than a threshold." << endl;
        cout << "If you are using prep file, this may take a while, please wait ..." << endl;
    }

    PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* recognized_residue = *it;
        string residue_name = recognized_residue->GetResidueName();
        if(residue_name.compare("HIS") != 0)
        {
            char chain_id = recognized_residue->GetResidueChainId();
            int sequence_number = recognized_residue->GetResidueSequenceNumber();
            char insertion_code = recognized_residue->GetResidueInsertionCode();
            char alternate_location = recognized_residue->GetResidueAlternateLocation();

            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            PdbFile::PdbAtomVector atoms_of_residue = *(residue_atom_map[key]);
            // Slow version
//            vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
            // Advanced version
            AtomNameMap dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFilesMap(residue_name, lib_files, prep_files);

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

void PdbPreprocessor::RemoveResiduesOfUnknownHeavyAtoms(PdbFile *pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms)
{
    vector<string> removed_keys = vector<string>();
    for(PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unknown_heavy_atoms.begin(); it != unknown_heavy_atoms.end(); it++)
    {
        PdbPreprocessorUnrecognizedHeavyAtom* unknown_heavy_atom = (*it);
        stringstream ss;
        ss << unknown_heavy_atom->GetResidueInsertionCode() << "_" << unknown_heavy_atom->GetResidueChainId() << "_" << unknown_heavy_atom->GetResidueSequenceNumber()
           << "_" << unknown_heavy_atom->GetResidueInsertionCode() << "_" << unknown_heavy_atom->GetResidueAlternateLocation();
        string residue_key = ss.str();
        if(distance(removed_keys.begin(), find(removed_keys.begin(), removed_keys.end(), residue_key)) < 0 ||
                distance(removed_keys.begin(), find(removed_keys.begin(), removed_keys.end(), residue_key)) >= (int)removed_keys.size())
        {
            PdbResidue* pdb_residue = new PdbResidue(unknown_heavy_atom->GetResidueName(),unknown_heavy_atom->GetResidueChainId(), unknown_heavy_atom->GetResidueSequenceNumber(),
                                                     unknown_heavy_atom->GetResidueInsertionCode(), unknown_heavy_atom->GetResidueAlternateLocation());
            pdb_file->DeleteResidue(pdb_residue);
            removed_keys.push_back(residue_key);
        }
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

vector<string> PdbPreprocessor::GetRemovedHydrogenNamesOfResidue(vector<string> pdb_atom_names_of_residue, AtomNameMap dataset_atom_names_of_residue)
{
    vector<string> removed_hydrogen_names_of_residue;
    for(vector<string>::iterator it = pdb_atom_names_of_residue.begin(); it != pdb_atom_names_of_residue.end(); it++)
    {
        string pdb_atom_name = (*it);
        if((pdb_atom_name.substr(0,1).compare("H") == 0 ||
            (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(dataset_atom_names_of_residue.find(pdb_atom_name) == dataset_atom_names_of_residue.end())
            {
                removed_hydrogen_names_of_residue.push_back(pdb_atom_name);
            }
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

PdbFileSpace::PdbFile::PdbAtomVector PdbPreprocessor::GetRemovedHydrogensOfResidue(PdbFile::PdbAtomVector pdb_atoms, AtomNameMap dataset_atom_names_of_residue)
{
    PdbFile::PdbAtomVector removed_hydrogens_of_residue;
    for(PdbFile::PdbAtomVector::iterator it = pdb_atoms.begin(); it != pdb_atoms.end(); it++)
    {
        PdbAtom* pdb_atom = *it;
        string pdb_atom_name = pdb_atom->GetAtomName();

        if((pdb_atom_name.substr(0,1).compare("H") == 0 ||
            (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(dataset_atom_names_of_residue.find(pdb_atom_name) == dataset_atom_names_of_residue.end())
                removed_hydrogens_of_residue.push_back(pdb_atom);
        }
    }
    return removed_hydrogens_of_residue;
}
void PdbPreprocessor::ExtractRemovedHydrogens(string pdb_file_path, vector<string> lib_files, vector<string> prep_files)
{
    try
    {
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        // Slow version
//        vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
        // Advanced version
        ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
        vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
        // Slow version
//        vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
        // Advanced version
        ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

        if(recognized_residues.size() > PdbResidueThreshold)
        {
            cout << "Number of recognized residues in this pdb file is more than a threshold." << endl;
            cout << "If you are using prep file, this may take a while, please wait ..." << endl;
        }

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
            // Slow version
//            vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
            // Advanced version
            AtomNameMap dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFilesMap(residue_name, lib_files, prep_files);
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
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractRemovedHydrogens(PdbFile* pdb_file, vector<string> lib_files, vector<string> prep_files)
{
    // Slow version
//    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    // Advanced version
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    // Slow version
//    vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    // Advanced version
    ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

    if(recognized_residues.size() > PdbResidueThreshold)
    {
        cout << "Number of recognized residues in this pdb file is more than a threshold." << endl;
        cout << "If you are using prep file, this may take a while, please wait ..." << endl;
    }

    PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* recognized_residue = *it;
        if(recognized_residue->GetResidueName().compare("HIS") != 0)
        {
            string residue_name = recognized_residue->GetResidueName();
            char chain_id = recognized_residue->GetResidueChainId();
            int sequence_number = recognized_residue->GetResidueSequenceNumber();
            char insertion_code = recognized_residue->GetResidueInsertionCode();
            char alternate_location = recognized_residue->GetResidueAlternateLocation();
            stringstream ss;
            ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            string key = ss.str();
            PdbFile::PdbAtomVector atoms_of_residue = *(residue_atom_map[key]);
            // Slow version
//            vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
            // Advanced version
            AtomNameMap dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFilesMap(residue_name, lib_files, prep_files);
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
    try
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
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractAminoAcidChains(PdbFile* pdb_file)
{
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
    try
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
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractGapsInAminoAcidChains(PdbFile* pdb_file)
{
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
void PdbPreprocessor::ExtractAlternateResidue(string pdb_file_path)
{
    try
    {
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
        {
            PdbFileSpace::PdbResidue* target_residue = *it;
            string target_residue_name = target_residue->GetResidueName();
            char target_chain_id = target_residue->GetResidueChainId();
            int target_sequence_number = target_residue->GetResidueSequenceNumber();
            char target_insertion_code = target_residue->GetResidueInsertionCode();
            char target_alternate_location = target_residue->GetResidueAlternateLocation();
            stringstream ss;
            ss << target_residue_name << "_" << target_chain_id << "_" << target_sequence_number << "_" << target_insertion_code;
            string target_key = ss.str();
            for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 =it+1 ; it1 != pdb_residues.end(); it1++)
            {
                PdbFileSpace::PdbResidue* residue = *it1;
                string residue_name = residue->GetResidueName();
                char chain_id = residue->GetResidueChainId();
                int sequence_number = residue->GetResidueSequenceNumber();
                char insertion_code = residue->GetResidueInsertionCode();
                char alternate_location = residue->GetResidueAlternateLocation();
                stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code;
                string key = sss.str();

                if(key.compare(target_key) == 0)
                {
                    if(target_alternate_location != alternate_location)
                    {
                        if (alternate_residue_map_.empty() || distance(alternate_residue_map_.begin() ,alternate_residue_map_.find(key)) < 0 ||
                                distance(alternate_residue_map_.begin() ,alternate_residue_map_.find(key)) >= (int)alternate_residue_map_.size())
                        {
                            vector<bool> selected = vector<bool>();
                            selected.push_back(true);
                            selected.push_back(false);
                            vector<char> alternate_locations = vector<char>();
                            alternate_locations.push_back(target_alternate_location);
                            alternate_locations.push_back(alternate_location);
                            alternate_residue_map_[target_key] = new PdbPreprocessorAlternateResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_locations, selected);
                        }
                        else
                        {
                            PdbPreprocessorAlternateResidue* alternate_residue = alternate_residue_map_[target_key];
                            vector<char> alternate_locations = alternate_residue->GetResidueAlternateLocation();
                            if(distance(alternate_locations.begin(), find(alternate_locations.begin(),alternate_locations.end(),alternate_location)) < 0 ||
                                    distance(alternate_locations.begin(), find(alternate_locations.begin(),alternate_locations.end(),alternate_location)) >= (int)alternate_locations.size())
                            {

                                alternate_locations.push_back(alternate_location);
                                alternate_residue->SetResidueAlternateLocation(alternate_locations);
                                vector<bool> selected_alternate_locations = alternate_residue->GetSelectedAlternateLocation();
                                selected_alternate_locations.push_back(false);
                                alternate_residue->SetSelectedAlternateLocation(selected_alternate_locations);
                            }
                        }
                    }
                }
            }
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::ExtractAlternateResidue(PdbFile* pdb_file)
{
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* target_residue = *it;
        string target_residue_name = target_residue->GetResidueName();
        char target_chain_id = target_residue->GetResidueChainId();
        int target_sequence_number = target_residue->GetResidueSequenceNumber();
        char target_insertion_code = target_residue->GetResidueInsertionCode();
        char target_alternate_location = target_residue->GetResidueAlternateLocation();
        stringstream ss;
        ss << target_residue_name << "_" << target_chain_id << "_" << target_sequence_number << "_" << target_insertion_code;
        string target_key = ss.str();
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 =it+1 ; it1 != pdb_residues.end(); it1++)
        {
            PdbFileSpace::PdbResidue* residue = *it1;
            string residue_name = residue->GetResidueName();
            char chain_id = residue->GetResidueChainId();
            int sequence_number = residue->GetResidueSequenceNumber();
            char insertion_code = residue->GetResidueInsertionCode();
            char alternate_location = residue->GetResidueAlternateLocation();
            stringstream sss;
            sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code;
            string key = sss.str();

            if(key.compare(target_key) == 0)
            {
                if(target_alternate_location != alternate_location)
                {
                    if (alternate_residue_map_.empty() || distance(alternate_residue_map_.begin() ,alternate_residue_map_.find(key)) < 0 ||
                            distance(alternate_residue_map_.begin() ,alternate_residue_map_.find(key)) >= (int)alternate_residue_map_.size())
                    {
                        vector<bool> selected = vector<bool>();
                        selected.push_back(true);
                        selected.push_back(false);
                        vector<char> alternate_locations = vector<char>();
                        alternate_locations.push_back(target_alternate_location);
                        alternate_locations.push_back(alternate_location);
                        alternate_residue_map_[target_key] = new PdbPreprocessorAlternateResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_locations, selected);
                    }
                    else
                    {
                        PdbPreprocessorAlternateResidue* alternate_residue = alternate_residue_map_[target_key];
                        vector<char> alternate_locations = alternate_residue->GetResidueAlternateLocation();
                        if(distance(alternate_locations.begin(), find(alternate_locations.begin(),alternate_locations.end(),alternate_location)) < 0 ||
                                distance(alternate_locations.begin(), find(alternate_locations.begin(),alternate_locations.end(),alternate_location)) >= (int)alternate_locations.size())
                        {

                            alternate_locations.push_back(alternate_location);
                            alternate_residue->SetResidueAlternateLocation(alternate_locations);
                            vector<bool> selected_alternate_locations = alternate_residue->GetSelectedAlternateLocation();
                            selected_alternate_locations.push_back(false);
                            alternate_residue->SetSelectedAlternateLocation(selected_alternate_locations);
                        }
                    }
                }
            }
        }
    }
}
void PdbPreprocessor::RemoveUnselectedAlternateResidues(PdbFile *pdb_file, PdbPreprocessorAlternateResidueMap alternate_residue_map)
{
    for(PdbPreprocessorAlternateResidueMap::iterator it = alternate_residue_map.begin(); it != alternate_residue_map.end(); it++)
    {
        PdbPreprocessorAlternateResidue* alternate_residue = (*it).second;
        vector<bool> selected_alternate_locations = alternate_residue->GetSelectedAlternateLocation();
        for(vector<bool>::iterator it1 = selected_alternate_locations.begin(); it1 != selected_alternate_locations.end(); it1++)
        {
            bool selected_alternate_location = (*it1);
            char alternate_location = alternate_residue->GetResidueAlternateLocation().at(distance(selected_alternate_locations.begin(), it1));
            if(!selected_alternate_location)
            {
                PdbResidue* pdb_residue = new PdbResidue(alternate_residue->GetResidueName() ,alternate_residue->GetResidueChainId(), alternate_residue->GetResidueSequenceNumber(),
                                                         alternate_residue->GetResidueInsertionCode(), alternate_location);
                pdb_file->DeleteResidue(pdb_residue);
            }
        }
    }

}

void PdbPreprocessor::DeleteAllToBeDeletedEntities(PdbFile *pdb_file)
{
    for(PdbPreprocessorToBeDeletedAtomVector::iterator it = to_be_deleted_atoms_.begin(); it != to_be_deleted_atoms_.end(); it++)
    {
        PdbAtom* atom = (*it);
        pdb_file->DeleteAtom(atom);

    }
    for(PdbPreprocessorToBeDeletedResidueVector::iterator it = to_be_deleted_residues_.begin(); it != to_be_deleted_residues_.end(); it++)
    {
        PdbResidue* residue = (*it);
        pdb_file->DeleteResidue(residue);

    }

}

void PdbPreprocessor::Preprocess(PdbFile* pdb_file, vector<string> lib_files_path, vector<string> prep_files_path)
{
    try
    {
        time_t t = time(0);
        cout << std::asctime(std::localtime(&t)) << "Start preprocessing ..." << endl;
        ExtractHISResidues(pdb_file);
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "HIS residues extraction: done" << endl;
        ExtractCYSResidues(pdb_file);
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "CYS residues extraction: done" << endl;
        ExtractAlternateResidue(pdb_file);
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "Alternate residues extraction: done" << endl;
        ExtractUnrecognizedResidues(pdb_file, lib_files_path, prep_files_path);
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "Unrecognized residues extraction: done" << endl;
        ExtractUnknownHeavyAtoms(pdb_file, lib_files_path, prep_files_path);
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "Unknown heavy atoms extraction: done" << endl;
        ExtractRemovedHydrogens(pdb_file, lib_files_path, prep_files_path);
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "Removed hydrogens extraction: done" << endl;
        ExtractAminoAcidChains(pdb_file);
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "Amino acid chains extraction: done" << endl;
        ExtractGapsInAminoAcidChains(pdb_file);
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "Gaps in amino acid chains extraction: done" << endl;
        t = time(0);
        cout << std::asctime(std::localtime(&t)) << "Preprocessing done" << endl;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

void PdbPreprocessor::ApplyPreprocessing(PdbFile *pdb_file, vector<string> lib_files_path)
{
    time_t t = time(0);
    cout << std::asctime(std::localtime(&t)) << "Start to apply changes ..." << endl;
    UpdateHISMapping(pdb_file,this->GetHistidineMappings());
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "HIS residues update: done" << endl;
    UpdateCYSResidues(pdb_file, this->GetDisulfideBonds());
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "CYS residues update: done" << endl;
    RemoveUnselectedAlternateResidues(pdb_file,this->GetAlternateResidueMap());
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "Unselected alternate residues removed: done" << endl;
    RemoveUnrecognizedResidues(pdb_file, this->GetUnrecognizedResidues());
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "Remove unrecognized residues: done" << endl;
    RemoveResiduesOfUnknownHeavyAtoms(pdb_file, this->GetUnrecognizedHeavyAtoms());
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "Unknown heavy atoms removed: done" << endl;
    RemoveRemovedHydrogens(pdb_file, this->GetReplacedHydrogens());
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "Removed hydrogens removed: done" << endl;
    UpdateAminoAcidChains(pdb_file,lib_files_path, this->GetChainTerminations());
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "Amino acid chains update: done" << endl;
    UpdateGapsInAminoAcidChains(pdb_file, lib_files_path, this->GetMissingResidues());
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "Gaps in amino acid chains update: done" << endl;
    t = time(0);
    cout << std::asctime(std::localtime(&t)) << "Applying changes done" << endl;
}


//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessor::Print(ostream &out)
{
    out << "================================== Disulfide Bonds =====================================" << endl;
    for(PdbPreprocessorDisulfideBondVector::iterator it = disulfide_bonds_.begin(); it != disulfide_bonds_.end(); it++)
    {
        PdbPreprocessorDisulfideBond* disulfide_bond = (*it);
        disulfide_bond->Print(out);
    }
    out << "======================================= Chains =========================================" << endl;
    for(PdbPreprocessorChainTerminationVector::iterator it = chain_terminations_.begin(); it != chain_terminations_.end(); it++)
    {
        PdbPreprocessorChainTermination* chain = (*it);
        chain->Print(out);
    }
    out << "===================================== HIS Residues =====================================" << endl;
    for(PdbPreprocessorHistidineMappingVector::iterator it = histidine_mappings_.begin(); it != histidine_mappings_.end(); it++)
    {
        PdbPreprocessorHistidineMapping* his_residue = (*it);
        his_residue->Print(out);
    }
    out << "======================================== Gaps ==========================================" << endl;
    for(PdbPreprocessorMissingResidueVector::iterator it = missing_residues_.begin(); it != missing_residues_.end(); it++)
    {
        PdbPreprocessorMissingResidue* gap = (*it);
        gap->Print(out);
    }
    out << "============================== Unrecognized Residues ===================================" << endl;
    for(PdbPreprocessorUnrecognizedResidueVector::iterator it = unrecognized_residues_.begin(); it != unrecognized_residues_.end(); it++)
    {
        PdbPreprocessorUnrecognizedResidue* unrecognized_residue = (*it);
        unrecognized_residue->Print(out);
    }
    out << "=============================== Unknown Heavy Atoms ====================================" << endl;
    for(PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unrecognized_heavy_atoms_.begin(); it != unrecognized_heavy_atoms_.end(); it++)
    {
        PdbPreprocessorUnrecognizedHeavyAtom* unknown_heavy_atom = (*it);
        unknown_heavy_atom->Print(out);
    }
    out << "================================ Removed Hydrogens =====================================" << endl;
    for(PdbPreprocessorReplacedHydrogenVector::iterator it = replaced_hydrogens_.begin(); it != replaced_hydrogens_.end(); it++)
    {
        PdbPreprocessorReplacedHydrogen* removed_hydrogen = (*it);
        removed_hydrogen->Print(out);
    }
    out << "================================ Alternate Residues ====================================" << endl;
    for(PdbPreprocessorAlternateResidueMap::iterator it = alternate_residue_map_.begin(); it != alternate_residue_map_.end(); it++)
    {
        PdbPreprocessorAlternateResidue* alternate_residue = (*it).second;
        alternate_residue->Print(out);
    }
}

