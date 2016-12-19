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
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
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
PdbPreprocessor::PdbPreprocessorResidueInfoMap PdbPreprocessor::GetResidueInfoMap(){
    return residue_info_map_;
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
vector<string> PdbPreprocessor::GetUnrecognizedResidueNames(PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, vector<string> dataset_residue_names)
{
    vector<string> unrecognized_residue_names;
    pair<string, string> residue_name_flag;
    string pdb_residue_name;
    string pdb_residue_position_flag;
    for(PdbFile::PdbPairVectorAtomNamePositionFlag::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        residue_name_flag = *it;
        pdb_residue_name = residue_name_flag.first;
        pdb_residue_position_flag = residue_name_flag.second;
        if(pdb_residue_name.compare("HIS") != 0 )
        {
            if(find(dataset_residue_names.begin(), dataset_residue_names.end(), pdb_residue_name) != dataset_residue_names.end())
            {
            }
            else
            {
                if(pdb_residue_position_flag.compare("S") == 0)
                {
                    stringstream ss;
                    ss << "N" << pdb_residue_name;
                    string n_terminal_pdb_residue_name = ss.str();
                    if(find(dataset_residue_names.begin(), dataset_residue_names.end(), n_terminal_pdb_residue_name) == dataset_residue_names.end())
                        unrecognized_residue_names.push_back(pdb_residue_name);
                }
                else if(pdb_residue_position_flag.compare("E") == 0)
                {
                    stringstream ss;
                    ss << "C" << pdb_residue_name;
                    string c_terminal_pdb_residue_name = ss.str();
                    if(find(dataset_residue_names.begin(), dataset_residue_names.end(), c_terminal_pdb_residue_name) == dataset_residue_names.end())
                        unrecognized_residue_names.push_back(pdb_residue_name);
                }
                else
                {
                    unrecognized_residue_names.push_back(pdb_residue_name);
                }
            }
        }
    }
    cout << "HIS residue(s) found" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "HIS residue(s) found" );
    return unrecognized_residue_names;
}

ResidueNameMap PdbPreprocessor::GetUnrecognizedResidueNamesMap(PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, ResidueNameMap dataset_residue_names)
{
    ResidueNameMap unrecognized_residue_names = ResidueNameMap();
    pair<string, string> residue_name_flag;
    string pdb_residue_name;
    string pdb_residue_position_flag;
    for(PdbFile::PdbPairVectorAtomNamePositionFlag::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
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
                    stringstream ss;
                    ss << "N" << pdb_residue_name;
                    string n_terminal_pdb_residue_name = ss.str();
                    if(dataset_residue_names.find(n_terminal_pdb_residue_name) == dataset_residue_names.end())
                        unrecognized_residue_names[pdb_residue_name] = pdb_residue_name;
                }
                else if(pdb_residue_position_flag.compare("E") == 0)
                {
                    stringstream ss;
                    ss << "C" << pdb_residue_name;
                    string c_terminal_pdb_residue_name = ss.str();
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
    cout << "HIS residue(s) found" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "HIS residue(s) found" );
    return unrecognized_residue_names;
}

vector<string> PdbPreprocessor::GetRecognizedResidueNames(PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, vector<string> dataset_residue_names)
{
    vector<string> recognized_residue_names;
    pair<string, string> residue_name_flag;
    string pdb_residue_name;
    string pdb_residue_position_flag;
    for(PdbFile::PdbPairVectorAtomNamePositionFlag::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
    {
        residue_name_flag = *it;
        pdb_residue_name = residue_name_flag.first;
        pdb_residue_position_flag = residue_name_flag.second;
        if(pdb_residue_name.compare("HIS") != 0)
        {
            if(find(dataset_residue_names.begin(), dataset_residue_names.end(), pdb_residue_name) != dataset_residue_names.end())
            {
                recognized_residue_names.push_back(pdb_residue_name);
            }
            else
            {
                if(pdb_residue_position_flag.compare("S") == 0)
                {
                    stringstream ss;
                    ss << "N" << pdb_residue_name;
                    string n_terminal_pdb_residue_name = ss.str();
                    if(find(dataset_residue_names.begin(), dataset_residue_names.end(), n_terminal_pdb_residue_name) != dataset_residue_names.end())
                        recognized_residue_names.push_back(pdb_residue_name);
                }
                else if(pdb_residue_position_flag.compare("E") == 0)
                {
                    stringstream ss;
                    ss << "C" << pdb_residue_name;
                    string c_terminal_pdb_residue_name = ss.str();
                    if(find(dataset_residue_names.begin(), dataset_residue_names.end(), c_terminal_pdb_residue_name) != dataset_residue_names.end())
                        recognized_residue_names.push_back(pdb_residue_name);
                }
            }
        }
        else
        {
            recognized_residue_names.push_back(pdb_residue_name);
        }
    }
    cout << "HIS residue(s) found" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "HIS residue(s) found" );
    return recognized_residue_names;
}

ResidueNameMap PdbPreprocessor::GetRecognizedResidueNamesMap(PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, ResidueNameMap dataset_residue_names)
{
    ResidueNameMap recognized_residue_names = ResidueNameMap();
    pair<string, string> residue_name_flag;
    string pdb_residue_name;
    string pdb_residue_position_flag;
    for(PdbFile::PdbPairVectorAtomNamePositionFlag::iterator it = pdb_residue_names.begin(); it != pdb_residue_names.end(); it++)
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
                    stringstream ss;
                    ss << "N" << pdb_residue_name;
                    string n_terminal_pdb_residue_name = ss.str();
                    if(dataset_residue_names.find(n_terminal_pdb_residue_name) != dataset_residue_names.end())
                        recognized_residue_names[pdb_residue_name] = pdb_residue_name;
                }
                else if(pdb_residue_position_flag.compare("E") == 0)
                {
                    stringstream ss;
                    ss << "C" << pdb_residue_name;
                    string c_terminal_pdb_residue_name = ss.str();
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
    cout << "HIS residue(s) found" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "HIS residue(s) found" );
    return recognized_residue_names;
}

PdbFileSpace::PdbFile::PdbResidueVector PdbPreprocessor::GetUnrecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, vector<string> unrecognized_residue_names)
{
    PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* pdb_residue = *it;
        string pdb_residue_name = pdb_residue->GetResidueName();
        if(find(unrecognized_residue_names.begin(), unrecognized_residue_names.end(), pdb_residue_name) == unrecognized_residue_names.end())
        {
            unrecognized_residues.push_back(pdb_residue);
        }
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
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* pdb_residue = *it;
        string pdb_residue_name = pdb_residue->GetResidueName();
        if(find(recognized_residue_names.begin(), recognized_residue_names.end(), pdb_residue_name) != recognized_residue_names.end())
            recognized_residues.push_back(pdb_residue);
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
    vector<string> residue_names;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        residue_names = lib_file->GetAllResidueNames();
        for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
        {
            string residue_name = (*it1);
            all_residue_names[residue_name] = (residue_name);
        }
    }
    return all_residue_names;
}

LibraryFileSpace::LibraryFile::ResidueMap PdbPreprocessor::GetAllResiduesFromMultipleLibFilesMap(vector<string> lib_files)
{
    LibraryFileSpace::LibraryFile::ResidueMap all_residues;
    LibraryFileSpace::LibraryFile::ResidueMap residues;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        residues = lib_file->GetResidues();
        for(LibraryFileSpace::LibraryFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            string lib_residue_name = (*it1).first;
            all_residues[lib_residue_name] = (*it1).second;
        }
    }
    return all_residues;
}
PrepFileSpace::PrepFile::ResidueMap PdbPreprocessor::GetAllResiduesFromMultiplePrepFilesMap(vector<string> prep_files)
{
    PrepFileSpace::PrepFile::ResidueMap all_residues;
    PrepFileSpace::PrepFile::ResidueMap residues;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        residues = prep_file->GetResidues();
        for(PrepFileSpace::PrepFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            string prep_residue_name = (*it1).first;
            all_residues[prep_residue_name] = (*it1).second;
        }
    }
    return all_residues;
}

vector<string> PdbPreprocessor::GetAllResidueNamesFromMultiplePrepFiles(vector<string> prep_files)
{
    vector<string> all_residue_names;
    vector<string> residue_names;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        residue_names = prep_file->GetAllResidueNames();
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
    vector<string> residue_names;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        residue_names = prep_file->GetAllResidueNames();
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

bool PdbPreprocessor::ExtractUnrecognizedResidues(string pdb_file_path, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    try
    {
        vector<string> lib_files = vector<string>();
        for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        // Slow version
//        vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
        // Advanced version
        ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = pdb_file->GetAllResidueNames();
        // Slow version
//        vector<string> unrecognized_residue_names = GetUnrecognizedResidueNames(pdb_residue_names, dataset_residue_names);
        // Advanced version
        ResidueNameMap unrecognized_residue_names = GetUnrecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);

        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues = GetUnrecognizedResidues(pdb_residues, unrecognized_residue_names);

        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues_atom_card = pdb_file->GetAllResiduesFromAtomCard();
        PdbPreprocessorChainIdResidueMap all_chain_map_residue;
        PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
        PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues_atom_card.begin(); it != pdb_residues_atom_card.end(); it++)
        {
            PdbFileSpace::PdbResidue* residue = *it;
            char chain_id = residue->GetResidueChainId();

            stringstream ss;
            ss << chain_id;
            all_chain_map_residue[ss.str()].push_back(residue);
        }

        // Slower version
    //        vector<string> amino_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
        // Advanced version
        ResidueNameMap amino_lib_residue_names_map = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);
        PdbPreprocessorChainIdResidueMap chain_map_residue;
        PdbFileSpace::PdbFile::PdbResidueVector residues;
        PdbFile::PdbResidueVector pdb_residues_in_chain;
        for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
        {
            int internal_amino_acid_chain_counter = 0;
            chain_map_residue = PdbPreprocessorChainIdResidueMap();
            string chain_id = (*it).first;
            residues = (*it).second;

            for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
            {
                PdbResidue* residue = (*it1);
                if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
                {
                    stringstream ss;
                    ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    stringstream sss(chain_id);
                    char c_id;
                    sss >> c_id;
                    chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                    chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
                }
                else
                {
                    stringstream ss;
                    ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    internal_amino_acid_chain_counter++;
                }
            }

            for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = unrecognized_residues.begin(); it1 != unrecognized_residues.end(); it1++)
            {
                PdbResidue* pdb_residue = (*it1);
                PdbPreprocessorUnrecognizedResidue* unrecognized_residue =
                        new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(), pdb_residue->GetResidueSequenceNumber(),
                                                               pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation(), false);
                if(chain_map_residue.size() <= 2)
                    unrecognized_residue->SetMiddleOfChain(false);
                else
                {
                    for(PdbPreprocessorChainIdResidueMap::iterator it2 = chain_map_residue.begin(); it2 != chain_map_residue.end(); it2++)
                    {
                        string key = (*it2).first;
                        pdb_residues_in_chain = (*it2).second;
                        int dist = distance(chain_map_residue.begin(), it2);
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
                unrecognized_residues_.push_back(unrecognized_residue);
            }
        }
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractUnrecognizedResidues(PdbFile* pdb_file, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    vector<string> lib_files = vector<string>();
    for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    // Slow version
//    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    // Advanced version
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
    PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = pdb_file->GetAllResidueNames();
    // Slow version
//    vector<string> unrecognized_residue_names = GetUnrecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    // Advanced version
    ResidueNameMap unrecognized_residue_names = GetUnrecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);

    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector unrecognized_residues = GetUnrecognizedResidues(pdb_residues, unrecognized_residue_names);

    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues_atom_card = pdb_file->GetAllResiduesFromAtomCard();
    PdbPreprocessorChainIdResidueMap all_chain_map_residue;
    PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
    PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues_atom_card.begin(); it != pdb_residues_atom_card.end(); it++)
    {
        PdbFileSpace::PdbResidue* residue = *it;
        char chain_id = residue->GetResidueChainId();

        stringstream ss;
        ss << chain_id;
        all_chain_map_residue[ss.str()].push_back(residue);
    }

    // Slower version
//        vector<string> amino_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
    // Advanced version
    ResidueNameMap amino_lib_residue_names_map = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);
    PdbPreprocessorChainIdResidueMap chain_map_residue;
    PdbFileSpace::PdbFile::PdbResidueVector residues;
    PdbFile::PdbResidueVector pdb_residues_in_chain;
    for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
    {
        int internal_amino_acid_chain_counter = 0;
        chain_map_residue = PdbPreprocessorChainIdResidueMap();
        string chain_id = (*it).first;
        residues = (*it).second;

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            PdbResidue* residue = (*it1);
            if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
            {
                stringstream ss;
                ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                chain_map_residue[ss.str()].push_back(residue);
                stringstream sss(chain_id);
                char c_id;
                sss >> c_id;
                chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
            }
            else
            {
                stringstream ss;
                ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                chain_map_residue[ss.str()].push_back(residue);
                internal_amino_acid_chain_counter++;
            }
        }

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = unrecognized_residues.begin(); it1 != unrecognized_residues.end(); it1++)
        {
            PdbResidue* pdb_residue = (*it1);
            PdbPreprocessorUnrecognizedResidue* unrecognized_residue =
                    new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(), pdb_residue->GetResidueSequenceNumber(),
                                                           pdb_residue->GetResidueInsertionCode(), pdb_residue->GetResidueAlternateLocation(), false);
            if(chain_map_residue.size() <= 2)
                unrecognized_residue->SetMiddleOfChain(false);
            else
            {
                for(PdbPreprocessorChainIdResidueMap::iterator it2 = chain_map_residue.begin(); it2 != chain_map_residue.end(); it2++)
                {
                    string key = (*it2).first;
                    pdb_residues_in_chain = (*it2).second;
                    int dist = distance(chain_map_residue.begin(), it2);
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
            unrecognized_residues_.push_back(unrecognized_residue);
        }
    }
    return true;
}
void PdbPreprocessor::RemoveUnrecognizedResidues(PdbFile *pdb_file, PdbPreprocessorUnrecognizedResidueVector unrecognized_residues)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorUnrecognizedResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
    {
        PdbPreprocessorUnrecognizedResidue* unrecognized_residue = (*it);
        PdbResidue* pdb_residue =
                new PdbResidue(unrecognized_residue->GetResidueName(), unrecognized_residue->GetResidueChainId(),
                               unrecognized_residue->GetResidueSequenceNumber(), unrecognized_residue->GetResidueInsertionCode(), unrecognized_residue->GetResidueAlternateLocation());
//        pdb_file->DeleteResidue(pdb_residue);
        to_be_deleted_residues_.push_back(pdb_residue);
    }
    DeleteAllToBeDeletedEntities(pdb_file);
}
void PdbPreprocessor::RemoveUnrecognizedResiduesWithTheGivenModelNumber(PdbFile *pdb_file, PdbPreprocessorUnrecognizedResidueVector unrecognized_residues, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorUnrecognizedResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
    {
        PdbPreprocessorUnrecognizedResidue* unrecognized_residue = (*it);
        PdbResidue* pdb_residue =
                new PdbResidue(unrecognized_residue->GetResidueName(), unrecognized_residue->GetResidueChainId(),
                               unrecognized_residue->GetResidueSequenceNumber(), unrecognized_residue->GetResidueInsertionCode(), unrecognized_residue->GetResidueAlternateLocation());
//        pdb_file->DeleteResidueWithTheGivenModelNumber(pdb_residue, model_number);
        to_be_deleted_residues_.push_back(pdb_residue);
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(pdb_file, model_number);
}
bool PdbPreprocessor::ExtractRecognizedResidues(string pdb_file_path, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    try
    {
        vector<string> lib_files = vector<string>();
        for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        // Slow version
//        vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
        // Advanced version
        ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = pdb_file->GetAllResidueNames();
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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractRecognizedResidues(PdbFile* pdb_file, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    vector<string> lib_files = vector<string>();
    for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    // Slow version
//    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    // Advanced version
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
    PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = pdb_file->GetAllResidueNames();
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
    return true;
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

double PdbPreprocessor::GetDistanceofCYS(PdbResidue *first_residue, PdbResidue *second_residue, PdbFile *pdb_file,
                                         PdbFile::PdbResidueAtomsMap residue_atom_map, int &first_sulfur_atom_serial_number,
                                         int &second_sulfur_atom_serial_number)
{
    double distance = 0.0;
    PdbAtom* first_residue_sulfur_atom = pdb_file->GetAtomOfResidueByName(first_residue, "SG", residue_atom_map);
    PdbAtom* second_residue_sulfur_atom = pdb_file->GetAtomOfResidueByName(second_residue, "SG", residue_atom_map);
    first_sulfur_atom_serial_number = first_residue_sulfur_atom->GetAtomSerialNumber();
    second_sulfur_atom_serial_number = second_residue_sulfur_atom->GetAtomSerialNumber();
    if(first_residue_sulfur_atom != NULL && second_residue_sulfur_atom != NULL)
        distance = first_residue_sulfur_atom->GetAtomOrthogonalCoordinate().Distance(second_residue_sulfur_atom->GetAtomOrthogonalCoordinate());
    return distance;
}
bool PdbPreprocessor::ExtractCYSResidues(string pdb_file_path)
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
                int first_sulfur_atom_serial_number;
                int second_sulfur_atom_serial_number;
                double distance = GetDistanceofCYS(first_residue, second_residue, pdb_file, residue_atom_map, first_sulfur_atom_serial_number,
                                                   second_sulfur_atom_serial_number);
                if (distance < dSulfurCutoff)
                {
                    PdbPreprocessorDisulfideBond* disulfide_bond =
                            new PdbPreprocessorDisulfideBond(first_residue->GetResidueChainId(), second_residue->GetResidueChainId(),
                                                             first_residue->GetResidueSequenceNumber(), second_residue->GetResidueSequenceNumber(),
                                                             distance, true, first_residue->GetResidueInsertionCode(), second_residue->GetResidueInsertionCode(),
                                                             first_residue->GetResidueAlternateLocation(), second_residue->GetResidueAlternateLocation(),
                                                             first_sulfur_atom_serial_number, second_sulfur_atom_serial_number);
                    disulfide_bonds_.push_back(disulfide_bond);
                }
            }
        }
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractCYSResidues(PdbFile* pdb_file)
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
            int first_sulfur_atom_serial_number;
            int second_sulfur_atom_serial_number;
            double distance = GetDistanceofCYS(first_residue, second_residue, pdb_file, residue_atom_map, first_sulfur_atom_serial_number,
                                               second_sulfur_atom_serial_number);
            if (distance < dSulfurCutoff)
            {
                PdbPreprocessorDisulfideBond* disulfide_bond =
                        new PdbPreprocessorDisulfideBond(first_residue->GetResidueChainId(), second_residue->GetResidueChainId(),
                                                         first_residue->GetResidueSequenceNumber(), second_residue->GetResidueSequenceNumber(),
                                                         distance, true, first_residue->GetResidueInsertionCode(), second_residue->GetResidueInsertionCode(),
                                                         first_residue->GetResidueAlternateLocation(), second_residue->GetResidueAlternateLocation(),
                                                         first_sulfur_atom_serial_number, second_sulfur_atom_serial_number);
                disulfide_bonds_.push_back(disulfide_bond);
            }
        }
    }
    return true;
}
void PdbPreprocessor::UpdateCYSResidues(PdbFile *pdb_file, PdbPreprocessorDisulfideBondVector disulfide_bonds)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    PdbFile::PdbResidueVector pdb_residues;
    for(PdbPreprocessorDisulfideBondVector::iterator it = disulfide_bonds.begin(); it != disulfide_bonds.end(); it++)
    {
        PdbPreprocessorDisulfideBond* disulfide_bond = (*it);
        if(disulfide_bond->GetIsBonded())
        {
            PdbAtom* pdb_atom_1 =
                    new PdbAtom(disulfide_bond->GetResidueChainId1(), "HG", "CYS",
                                disulfide_bond->GetResidueSequenceNumber1(), disulfide_bond->GetResidueInsertionCode1(), disulfide_bond->GetResidueAlternateLocation1());
//            pdb_file->DeleteAtom(pdb_atom_1);
            to_be_deleted_atoms_.push_back(pdb_atom_1);
            PdbAtom* pdb_atom_2 =
                    new PdbAtom(disulfide_bond->GetResidueChainId2(), "HG", "CYS",
                                disulfide_bond->GetResidueSequenceNumber2(), disulfide_bond->GetResidueInsertionCode2(), disulfide_bond->GetResidueAlternateLocation2());
//            pdb_file->DeleteAtom(pdb_atom_2);
            to_be_deleted_atoms_.push_back(pdb_atom_2);
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
            pdb_residues = pdb_file->GetAllResidues();
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
            PdbConnectCard* pdb_connect_card;
            if(pdb_file->GetConnectivities() == NULL)
            {
                pdb_connect_card = new PdbConnectCard();
            }
            else
            {
                pdb_connect_card = pdb_file->GetConnectivities();
            }
            PdbConnectCard::BondedAtomsSerialNumbersMap bonded_atoms_serial_numbers_map = pdb_connect_card->GetBondedAtomsSerialNumbers();
            if(bonded_atoms_serial_numbers_map.find(disulfide_bond->GetSulfurAtomSerialNumber1()) != bonded_atoms_serial_numbers_map.end())
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()].push_back(disulfide_bond->GetSulfurAtomSerialNumber2());
            else
            {
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()] = vector<int>();
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()].push_back(disulfide_bond->GetSulfurAtomSerialNumber2());
            }
            if(bonded_atoms_serial_numbers_map.find(disulfide_bond->GetSulfurAtomSerialNumber2()) != bonded_atoms_serial_numbers_map.end())
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()].push_back(disulfide_bond->GetSulfurAtomSerialNumber1());
            else
            {
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()] = vector<int>();
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()].push_back(disulfide_bond->GetSulfurAtomSerialNumber1());
            }
            pdb_connect_card->SetBondedAtomsSerialNumbers(bonded_atoms_serial_numbers_map);
            pdb_file->SetConnectivities(pdb_connect_card);
        }
    }
    DeleteAllToBeDeletedEntities(pdb_file);
}

void PdbPreprocessor::UpdateCYSResiduesWithTheGivenModelNumber(PdbFile *pdb_file, PdbPreprocessorDisulfideBondVector disulfide_bonds, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    PdbFile::PdbResidueVector pdb_residues;
    for(PdbPreprocessorDisulfideBondVector::iterator it = disulfide_bonds.begin(); it != disulfide_bonds.end(); it++)
    {
        PdbPreprocessorDisulfideBond* disulfide_bond = (*it);
        if(disulfide_bond->GetIsBonded())
        {
            PdbAtom* pdb_atom_1 =
                    new PdbAtom(disulfide_bond->GetResidueChainId1(), "HG", "CYS",
                                disulfide_bond->GetResidueSequenceNumber1(), disulfide_bond->GetResidueInsertionCode1(), disulfide_bond->GetResidueAlternateLocation1());
//            pdb_file->DeleteAtomWithTheGivenModelNumber(pdb_atom_1, model_number);
            to_be_deleted_atoms_.push_back(pdb_atom_1);
            PdbAtom* pdb_atom_2 =
                    new PdbAtom(disulfide_bond->GetResidueChainId2(), "HG", "CYS",
                                disulfide_bond->GetResidueSequenceNumber2(), disulfide_bond->GetResidueInsertionCode2(), disulfide_bond->GetResidueAlternateLocation2());
//            pdb_file->DeleteAtomWithTheGivenModelNumber(pdb_atom_2, model_number);
            to_be_deleted_atoms_.push_back(pdb_atom_2);
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
            pdb_residues = pdb_file->GetAllResidues();
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
                        pdb_file->UpdateResidueNameWithTheGivenModelNumber(pdb_residue, "CYX", model_number);
                    }
                }
            }
            PdbConnectCard* pdb_connect_card;
            if(pdb_file->GetConnectivities() == NULL)
            {
                pdb_connect_card = new PdbConnectCard();
            }
            else
            {
                pdb_connect_card = pdb_file->GetConnectivities();
            }
            PdbConnectCard::BondedAtomsSerialNumbersMap bonded_atoms_serial_numbers_map = pdb_connect_card->GetBondedAtomsSerialNumbers();
            if(bonded_atoms_serial_numbers_map.find(disulfide_bond->GetSulfurAtomSerialNumber1()) != bonded_atoms_serial_numbers_map.end())
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()].push_back(disulfide_bond->GetSulfurAtomSerialNumber2());
            else
            {
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()] = vector<int>();
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber1()].push_back(disulfide_bond->GetSulfurAtomSerialNumber2());
            }
            if(bonded_atoms_serial_numbers_map.find(disulfide_bond->GetSulfurAtomSerialNumber2()) != bonded_atoms_serial_numbers_map.end())
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()].push_back(disulfide_bond->GetSulfurAtomSerialNumber1());
            else
            {
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()] = vector<int>();
                bonded_atoms_serial_numbers_map[disulfide_bond->GetSulfurAtomSerialNumber2()].push_back(disulfide_bond->GetSulfurAtomSerialNumber1());
            }

            pdb_connect_card->SetBondedAtomsSerialNumbers(bonded_atoms_serial_numbers_map);
            pdb_file->SetConnectivities(pdb_connect_card);
        }
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(pdb_file, model_number);
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
bool PdbPreprocessor::ExtractHISResidues(string pdb_file_path)
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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractHISResidues(PdbFile* pdb_file)
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
    return true;
}
void PdbPreprocessor::UpdateHISMapping(PdbFile *pdb_file, PdbPreprocessor::PdbPreprocessorHistidineMappingVector histidine_mappings)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
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
//                            pdb_file->DeleteAtom(HE2);
                            to_be_deleted_atoms_.push_back(HE2);
                        }
                    }
                    // HID residue
                    if(HE2 == NULL && HD1 != NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == HIE)
                        {
                            // Delete HD1
//                            pdb_file->DeleteAtom(HD1);
                            to_be_deleted_atoms_.push_back(HD1);
                        }

                    }
                    // HIP residue
                    if(HE2 != NULL && HD1 != NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == HIE)
                        {
                            // Delete HD1
                            pdb_file->DeleteAtom(HD1);
                            to_be_deleted_atoms_.push_back(HD1);
                        }
                        if(histidine_mapping->GetSelectedMapping() == HID)
                        {
                            // Delete HE2
//                            pdb_file->DeleteAtom(HE2);
                            to_be_deleted_atoms_.push_back(HE2);
                        }
                    }
                    pdb_file->UpdateResidueName(pdb_residue, histidine_mapping->GetStringFormatOfSelectedMapping());
                }
            }
        }
    }
    DeleteAllToBeDeletedEntities(pdb_file);
}

void PdbPreprocessor::UpdateHISMappingWithTheGivenNumber(PdbFile *pdb_file, PdbPreprocessorHistidineMappingVector histidine_mappings, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
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
//                            pdb_file->DeleteAtomWithTheGivenModelNumber(HE2, model_number);
                            to_be_deleted_atoms_.push_back(HE2);
                        }
                    }
                    // HID residue
                    if(HE2 == NULL && HD1 != NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == HIE)
                        {
                            // Delete HD1
//                            pdb_file->DeleteAtomWithTheGivenModelNumber(HD1, model_number);
                            to_be_deleted_atoms_.push_back(HD1);
                        }

                    }
                    // HIP residue
                    if(HE2 != NULL && HD1 != NULL)
                    {
                        if(histidine_mapping->GetSelectedMapping() == HIE)
                        {
                            // Delete HD1
//                            pdb_file->DeleteAtomWithTheGivenModelNumber(HD1, model_number);
                            to_be_deleted_atoms_.push_back(HD1);
                        }
                        if(histidine_mapping->GetSelectedMapping() == HID)
                        {
                            // Delete HE2
//                            pdb_file->DeleteAtomWithTheGivenModelNumber(HE2, model_number);
                            to_be_deleted_atoms_.push_back(HE2);
                        }
                    }
                    pdb_file->UpdateResidueNameWithTheGivenModelNumber(pdb_residue, histidine_mapping->GetStringFormatOfSelectedMapping(), model_number);
                }
            }
        }
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(pdb_file, model_number);
}

vector<string> PdbPreprocessor::GetUnknownHeavyAtomNamesOfResidue(vector<string> pdb_atom_names_of_residue, vector<string> dataset_atom_names_of_residue)
{
    vector<string> unknown_heavy_atom_names_of_residue;
    for(vector<string>::iterator it = pdb_atom_names_of_residue.begin(); it != pdb_atom_names_of_residue.end(); it++)
    {
        string pdb_atom_name = (*it);
        if(!(pdb_atom_name.substr(0,1).compare("H") == 0 ||
             (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(find(dataset_atom_names_of_residue.begin(), dataset_atom_names_of_residue.end(), pdb_atom_name) == dataset_atom_names_of_residue.end())
                unknown_heavy_atom_names_of_residue.push_back(pdb_atom_name);
        }
    }
    return unknown_heavy_atom_names_of_residue;
}

vector<string> PdbPreprocessor::GetAllAtomNamesOfResidueFromMultipleLibFiles(string residue_name, vector<string> lib_files)
{
    vector<string> all_atom_names_of_residue;
    vector<string> residue_names;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        residue_names = lib_file->GetAllResidueNames();
        if(find(residue_names.begin(), residue_names.end(), residue_name) != residue_names.end())
        {

            all_atom_names_of_residue = lib_file->GetAllAtomNamesOfResidue(residue_name);
            break;
        }
    }
    return all_atom_names_of_residue;
}

ResidueNameAtomNamesMap PdbPreprocessor::GetAllAtomNamesOfResidueNamesFromMultipleLibFiles(vector<string> lib_files)
{
    ResidueNameAtomNamesMap residue_atom_map = ResidueNameAtomNamesMap();
    LibraryFileSpace::LibraryFile::ResidueMap residues;
    for(vector<string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
    {
        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
        residues = lib_file->GetResidues();
        for(LibraryFileSpace::LibraryFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            string residue_name = (*it1).first;
            LibraryFileSpace::LibraryFileResidue* residue = (*it1).second;
            LibraryFileSpace::LibraryFileResidue::AtomMap atoms = residue->GetAtoms();
            for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                LibraryFileSpace::LibraryFileAtom* atom = (*it2).second;
                string atom_name = atom->GetName();
                residue_atom_map[residue_name].push_back(atom_name);
            }
        }
    }
    return residue_atom_map;
}

vector<string> PdbPreprocessor::GetAllAtomNamesOfResidueFromMultiplePrepFiles(string residue_name, vector<string> prep_files)
{
    vector<string> all_atom_names_of_residue;
    vector<string> residue_names;
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        residue_names = prep_file->GetAllResidueNames();
        if(find(residue_names.begin(), residue_names.end(), residue_name) != residue_names.end())
        {
            all_atom_names_of_residue = prep_file->GetAllAtomNamesOfResidue(residue_name);
            break;
        }
    }
    return all_atom_names_of_residue;
}

ResidueNameAtomNamesMap PdbPreprocessor::GetAllAtomNamesOfResidueNamesFromMultiplePrepFiles(vector<string> prep_files)
{
    ResidueNameAtomNamesMap residue_atom_map = ResidueNameAtomNamesMap();
    for(vector<string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
    {
        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
        PrepFileSpace::PrepFile::ResidueMap residues = prep_file->GetResidues();
        for(PrepFileSpace::PrepFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            string residue_name = (*it1).first;
            PrepFileSpace::PrepFileResidue* residue = (*it1).second;
            PrepFileSpace::PrepFileResidue::PrepFileAtomVector atoms = residue->GetAtoms();
            for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                PrepFileSpace::PrepFileAtom* atom = (*it2);
                string atom_name = atom->GetName();
                residue_atom_map[residue_name].push_back(atom_name);
            }
        }
    }
    return residue_atom_map;
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

ResidueNameAtomNamesMap PdbPreprocessor::GetAllAtomNamesOfResidueNamesFromDatasetFiles(vector<string> lib_files, vector<string> prep_files)
{
    ResidueNameAtomNamesMap all_residue_atom_map = ResidueNameAtomNamesMap();
    ResidueNameAtomNamesMap all_residue_atom_map_from_lib = GetAllAtomNamesOfResidueNamesFromMultipleLibFiles(lib_files);
    ResidueNameAtomNamesMap all_residue_atom_map_from_prep = GetAllAtomNamesOfResidueNamesFromMultiplePrepFiles(prep_files);
    for(ResidueNameAtomNamesMap::iterator it = all_residue_atom_map_from_lib.begin(); it != all_residue_atom_map_from_lib.end(); it++)
    {
        string residue_name = (*it).first;
        vector<string> atom_names = (*it).second;
        all_residue_atom_map[residue_name] = atom_names;
    }
    for(ResidueNameAtomNamesMap::iterator it = all_residue_atom_map_from_prep.begin(); it != all_residue_atom_map_from_prep.end(); it++)
    {
        string residue_name = (*it).first;
        vector<string> atom_names = (*it).second;
        all_residue_atom_map[residue_name] = atom_names;
    }
    return all_residue_atom_map;
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

bool PdbPreprocessor::ExtractUnknownHeavyAtoms(string pdb_file_path, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    try
    {
        vector<string> lib_files = vector<string>();
        for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        // Slow version
//         vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
        // Advanced version
        ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
        PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = pdb_file->GetAllResidueNames();
        // Slow version
//         vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
        // Advanced version
        ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);


        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

        PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
        ResidueNameAtomNamesMap dataset_residue_atom_map = GetAllAtomNamesOfResidueNamesFromDatasetFiles(lib_files, prep_files);
        vector<string> dataset_atom_names_of_residue = vector<string>();
        vector<string> dataset_atom_names_of_tail_residue = vector<string>();
        vector<string> dataset_atom_names_of_head_residue = vector<string>();
        PdbFile::PdbAtomVector unknown_heavy_atoms = PdbFile::PdbAtomVector();
        PdbFile::PdbAtomVector atoms_of_residue = PdbFile::PdbAtomVector();
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
                atoms_of_residue = *(residue_atom_map[key]);
                // Slow version
//                vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
                // Advanced version
                dataset_atom_names_of_residue = dataset_residue_atom_map[residue_name];
                pair<string, string> residue_sflag_pair = make_pair(residue_name, "S");
                pair<string, string> residue_eflag_pair = make_pair(residue_name, "E");
                if(find(pdb_residue_names.begin(), pdb_residue_names.end(), residue_sflag_pair) != pdb_residue_names.end())
                {
                    stringstream ss1;
                    ss1 << "N" << residue_name;
                    dataset_atom_names_of_head_residue = dataset_residue_atom_map[ss1.str()];
                    for(vector<string>::iterator it1 = dataset_atom_names_of_head_residue.begin(); it1 != dataset_atom_names_of_head_residue.end(); it1++)
                    {
                        string head_residue = (*it1);
                        dataset_atom_names_of_residue.push_back(head_residue);
                    }
                }
                else if(find(pdb_residue_names.begin(), pdb_residue_names.end(), residue_eflag_pair) != pdb_residue_names.end())
                {
                    stringstream ss1;
                    ss1 << "C" << residue_name;
                    dataset_atom_names_of_tail_residue = dataset_residue_atom_map[ss1.str()];
                    for(vector<string>::iterator it2 = dataset_atom_names_of_tail_residue.begin(); it2 != dataset_atom_names_of_tail_residue.end(); it2++)
                    {
                        string head_residue = (*it2);
                        dataset_atom_names_of_residue.push_back(head_residue);
                    }
                }
                unknown_heavy_atoms = GetUnknownHeavyAtomsOfResidue(atoms_of_residue, dataset_atom_names_of_residue);

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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractUnknownHeavyAtoms(PdbFile* pdb_file, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    vector<string> lib_files = vector<string>();
    for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    // Slow version
//    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    // Advanced version
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
    PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = pdb_file->GetAllResidueNames();
    // Slow version
//    vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    // Advanced version
    ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);


    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

    PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
    ResidueNameAtomNamesMap dataset_residue_atom_map = GetAllAtomNamesOfResidueNamesFromDatasetFiles(lib_files, prep_files);
    vector<string> dataset_atom_names_of_residue = vector<string>();
    vector<string> dataset_atom_names_of_head_residue = vector<string>();
    vector<string> dataset_atom_names_of_tail_residue = vector<string>();
    PdbFile::PdbAtomVector unknown_heavy_atoms = PdbFile::PdbAtomVector();
    PdbFile::PdbAtomVector atoms_of_residue = PdbFile::PdbAtomVector();
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
            atoms_of_residue = *(residue_atom_map[key]);
            // Slow version
//            vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
            // Advanced version
            dataset_atom_names_of_residue = dataset_residue_atom_map[residue_name];

            pair<string, string> residue_sflag_pair = make_pair(residue_name, "S");
            pair<string, string> residue_eflag_pair = make_pair(residue_name, "E");
            if(find(pdb_residue_names.begin(), pdb_residue_names.end(), residue_sflag_pair) != pdb_residue_names.end())
            {
                stringstream ss1;
                ss1 << "N" << residue_name;
                dataset_atom_names_of_head_residue = dataset_residue_atom_map[ss1.str()];
                for(vector<string>::iterator it1 = dataset_atom_names_of_head_residue.begin(); it1 != dataset_atom_names_of_head_residue.end(); it1++)
                {
                    string head_residue = (*it1);
                    dataset_atom_names_of_residue.push_back(head_residue);
                }
            }
            else if(find(pdb_residue_names.begin(), pdb_residue_names.end(), residue_eflag_pair) != pdb_residue_names.end())
            {
                stringstream ss1;
                ss1 << "C" << residue_name;
                dataset_atom_names_of_tail_residue = dataset_residue_atom_map[ss1.str()];
                for(vector<string>::iterator it2 = dataset_atom_names_of_tail_residue.begin(); it2 != dataset_atom_names_of_tail_residue.end(); it2++)
                {
                    string head_residue = (*it2);
                    dataset_atom_names_of_residue.push_back(head_residue);
                }
            }

            unknown_heavy_atoms = GetUnknownHeavyAtomsOfResidue(atoms_of_residue, dataset_atom_names_of_residue);

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
    return true;
}
void PdbPreprocessor::RemoveUnknownHeavyAtoms(PdbFile *pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unknown_heavy_atoms.begin(); it != unknown_heavy_atoms.end(); it++)
    {
        PdbPreprocessorUnrecognizedHeavyAtom* unknown_heavy_atom = (*it);
        PdbAtom* pdb_atom =
                new PdbAtom(unknown_heavy_atom->GetResidueChainId(),unknown_heavy_atom->GetAtomName(),
                            unknown_heavy_atom->GetResidueName(), unknown_heavy_atom->GetResidueSequenceNumber(), unknown_heavy_atom->GetResidueInsertionCode(), unknown_heavy_atom->GetResidueAlternateLocation());
//        pdb_file->DeleteAtom(pdb_atom);
        to_be_deleted_atoms_.push_back(pdb_atom);
    }
    DeleteAllToBeDeletedEntities(pdb_file);
}
void PdbPreprocessor::RemoveUnknownHeavyAtomsWithTheGivenModelNumber(PdbFile *pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms,int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unknown_heavy_atoms.begin(); it != unknown_heavy_atoms.end(); it++)
    {
        PdbPreprocessorUnrecognizedHeavyAtom* unknown_heavy_atom = (*it);
        PdbAtom* pdb_atom =
                new PdbAtom(unknown_heavy_atom->GetResidueChainId(),unknown_heavy_atom->GetAtomName(),
                            unknown_heavy_atom->GetResidueName(), unknown_heavy_atom->GetResidueSequenceNumber(), unknown_heavy_atom->GetResidueInsertionCode(), unknown_heavy_atom->GetResidueAlternateLocation());
//        pdb_file->DeleteAtomWithTheGivenModelNumber(pdb_atom, model_number);
        to_be_deleted_atoms_.push_back(pdb_atom);
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(pdb_file, model_number);
}

void PdbPreprocessor::RemoveResiduesOfUnknownHeavyAtoms(PdbFile *pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
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
//            pdb_file->DeleteResidue(pdb_residue);
            to_be_deleted_residues_.push_back(pdb_residue);
            removed_keys.push_back(residue_key);
        }
    }
    DeleteAllToBeDeletedEntities(pdb_file);
}

void PdbPreprocessor::RemoveResiduesOfUnknownHeavyAtomsWithTheGivenModelNumber(PdbFile *pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
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
//            pdb_file->DeleteResidueWithTheGivenModelNumber(pdb_residue , model_number);
            to_be_deleted_residues_.push_back(pdb_residue);
            removed_keys.push_back(residue_key);
        }
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(pdb_file, model_number);
}

vector<string> PdbPreprocessor::GetRemovedHydrogenNamesOfResidue(vector<string> pdb_atom_names_of_residue, vector<string> dataset_atom_names_of_residue)
{
    vector<string> removed_hydrogen_names_of_residue;
    for(vector<string>::iterator it = pdb_atom_names_of_residue.begin(); it != pdb_atom_names_of_residue.end(); it++)
    {
        string pdb_atom_name = (*it);
        if((pdb_atom_name.substr(0,1).compare("H") == 0 ||
            (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(find(dataset_atom_names_of_residue.begin(), dataset_atom_names_of_residue.end(), pdb_atom_name) == dataset_atom_names_of_residue.end())
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
        if((pdb_atom_name.substr(0,1).compare("H") == 0 ||
            (pdb_atom_name.substr(1,1).compare("H") == 0 && isdigit(ConvertString<char>(pdb_atom_name.substr(0,1))))))
        {
            if(find(dataset_atom_names_of_residue.begin(), dataset_atom_names_of_residue.end(), pdb_atom_name) == dataset_atom_names_of_residue.end())
                removed_hydrogens_of_residue.push_back(pdb_atom);
        }
    }
    return removed_hydrogens_of_residue;
}

bool PdbPreprocessor::ExtractRemovedHydrogens(string pdb_file_path, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    try
    {
        vector<string> lib_files = vector<string>();
        for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        // Slow version
//        vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
        // Advanced version
        ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
        PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = pdb_file->GetAllResidueNames();
        // Slow version
//        vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
        // Advanced version
        ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

        PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
        PdbFile::PdbAtomVector atoms_of_residue;
        ResidueNameAtomNamesMap dataset_residue_atom_map = GetAllAtomNamesOfResidueNamesFromDatasetFiles(lib_files, prep_files);
        vector<string> dataset_atom_names_of_residue = vector<string>();
        PdbFile::PdbAtomVector removed_hydrogens;
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
            atoms_of_residue = *(residue_atom_map[key]);
            // Slow version
//            vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
            // Advanced version
            dataset_atom_names_of_residue = dataset_residue_atom_map[residue_name];
            removed_hydrogens = GetRemovedHydrogensOfResidue(atoms_of_residue, dataset_atom_names_of_residue);
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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractRemovedHydrogens(PdbFile* pdb_file, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    vector<string> lib_files = vector<string>();
    for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    // Slow version
//    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    // Advanced version
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromDatasetFilesMap(lib_files, prep_files);
    PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = pdb_file->GetAllResidueNames();
    // Slow version
//    vector<string> recognized_residue_names = GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    // Advanced version
    ResidueNameMap recognized_residue_names = GetRecognizedResidueNamesMap(pdb_residue_names, dataset_residue_names);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
    PdbFileSpace::PdbFile::PdbResidueVector recognized_residues = GetRecognizedResidues(pdb_residues, recognized_residue_names);

    PdbFile::PdbResidueAtomsMap residue_atom_map = pdb_file->GetAllAtomsOfResidues();
    PdbFile::PdbAtomVector atoms_of_residue;
    ResidueNameAtomNamesMap dataset_residue_atom_map = GetAllAtomNamesOfResidueNamesFromDatasetFiles(lib_files, prep_files);
    vector<string> dataset_atom_names_of_residue = vector<string>();
    PdbFile::PdbAtomVector removed_hydrogens;
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
            atoms_of_residue = *(residue_atom_map[key]);
            // Slow version
//            vector<string> dataset_atom_names_of_residue = GetAllAtomNamesOfResidueFromDatasetFiles(residue_name, lib_files, prep_files);
            // Advanced version
            dataset_atom_names_of_residue = dataset_residue_atom_map[residue_name];
            removed_hydrogens = GetRemovedHydrogensOfResidue(atoms_of_residue, dataset_atom_names_of_residue);
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
    return true;
}
void PdbPreprocessor::RemoveRemovedHydrogens(PdbFile *pdb_file, PdbPreprocessorReplacedHydrogenVector replaced_hydrogens)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorReplacedHydrogenVector::iterator it = replaced_hydrogens.begin(); it != replaced_hydrogens.end(); it++)
    {
        PdbPreprocessorReplacedHydrogen* replaced_hydrogen = (*it);
        PdbAtom* pdb_atom =
                new PdbAtom(replaced_hydrogen->GetResidueChainId(),replaced_hydrogen->GetAtomName(),
                            replaced_hydrogen->GetResidueName(), replaced_hydrogen->GetResidueSequenceNumber(), replaced_hydrogen->GetResidueInsertionCode(), replaced_hydrogen->GetResidueAlternateLocation());
//        pdb_file->DeleteAtom(pdb_atom);
        to_be_deleted_atoms_.push_back(pdb_atom);
    }
    DeleteAllToBeDeletedEntities(pdb_file);
}
void PdbPreprocessor::RemoveRemovedHydrogensWithTheGivenModelNumber(PdbFile *pdb_file, PdbPreprocessorReplacedHydrogenVector replaced_hydrogens, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
    for(PdbPreprocessorReplacedHydrogenVector::iterator it = replaced_hydrogens.begin(); it != replaced_hydrogens.end(); it++)
    {
        PdbPreprocessorReplacedHydrogen* replaced_hydrogen = (*it);
        PdbAtom* pdb_atom =
                new PdbAtom(replaced_hydrogen->GetResidueChainId(),replaced_hydrogen->GetAtomName(),
                            replaced_hydrogen->GetResidueName(), replaced_hydrogen->GetResidueSequenceNumber(), replaced_hydrogen->GetResidueInsertionCode(), replaced_hydrogen->GetResidueAlternateLocation());
//        pdb_file->DeleteAtomWithTheGivenModelNumber(pdb_atom, model_number);
        to_be_deleted_atoms_.push_back(pdb_atom);
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(pdb_file, model_number);
}
bool PdbPreprocessor::ExtractAminoAcidChains(string pdb_file_path)
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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractAminoAcidChains(string pdb_file_path, vector<string> amino_lib_files)
{
    try
    {
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();

        PdbPreprocessorChainIdResidueMap all_chain_map_residue;
        PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
        PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
        {
            PdbFileSpace::PdbResidue* residue = *it;
            char chain_id = residue->GetResidueChainId();

            stringstream ss;
            ss << chain_id;
            all_chain_map_residue[ss.str()].push_back(residue);
        }

        // Slower version
//        vector<string> amino_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
        // Advanced version
        ResidueNameMap amino_lib_residue_names_map = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);
        amino_lib_residue_names_map["HIS"] = "HIS";

        for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
        {
            int internal_amino_acid_chain_counter = 0;
            PdbPreprocessorChainIdResidueMap chain_map_residue = PdbPreprocessorChainIdResidueMap();
            string chain_id = (*it).first;
            PdbFileSpace::PdbFile::PdbResidueVector residues = (*it).second;

            for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
            {
                PdbResidue* residue = (*it1);
                if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
                {
                    stringstream ss;
                    ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    stringstream sss(chain_id);
                    char c_id;
                    sss >> c_id;
                    chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                    chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
                }
                else
                {
                    stringstream ss;
                    ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    internal_amino_acid_chain_counter++;
                }
            }
            if(chain_map_residue.size() > 2)
            {
                cout << "There is an undefined protein in the middle of the chain" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is an undefined protein in the middle of the chain" );
                cout << "Pdb file is not processible at this time" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "Pdb file is not processible at this time" );

                return false;
            }
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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractAminoAcidChains(PdbFile* pdb_file)
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
    return true;
}
bool PdbPreprocessor::ExtractAminoAcidChains(PdbFile* pdb_file, vector<string> amino_lib_files)
{
    try
    {
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();

        PdbPreprocessorChainIdResidueMap all_chain_map_residue;
        PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
        PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
        {
            PdbFileSpace::PdbResidue* residue = *it;
            char chain_id = residue->GetResidueChainId();

            stringstream ss;
            ss << chain_id;
            all_chain_map_residue[ss.str()].push_back(residue);
        }

        // Slower version
//        vector<string> amino_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
        // Advanced version
        ResidueNameMap amino_lib_residue_names_map = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);
        amino_lib_residue_names_map["HIS"] = "HIS";

        for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
        {
            int internal_amino_acid_chain_counter = 0;
            PdbPreprocessorChainIdResidueMap chain_map_residue = PdbPreprocessorChainIdResidueMap();
            string chain_id = (*it).first;
            PdbFileSpace::PdbFile::PdbResidueVector residues = (*it).second;

            for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
            {
                PdbResidue* residue = (*it1);
                if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
                {
                    stringstream ss;
                    ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    stringstream sss(chain_id);
                    char c_id;
                    sss >> c_id;
                    chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                    chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
                }
                else
                {
                    stringstream ss;
                    ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    internal_amino_acid_chain_counter++;
                }
            }
            if(chain_map_residue.size() > 2)
            {
                cout << "There is an undefined protein in the middle of the chain" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is an undefined protein in the middle of the chain" );
                cout << "Pdb file is not processible at this time" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "Pdb file is not processible at this time" );


                return false;
            }
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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
void PdbPreprocessor::UpdateAminoAcidChains(PdbFile *pdb_file, vector<string> amino_lib_files, vector<string> glycam_lib_files,
                                            vector<string> prep_files, PdbPreprocessorChainTerminationVector chain_terminations)
{
    // Before non-amino-acid residue
    for(PdbPreprocessor::PdbPreprocessorChainTerminationVector::iterator it1 = chain_terminations.begin(); it1 != chain_terminations.end(); it1++)
    {
        PdbPreprocessorChainTermination* chain = (*it1);
        pdb_file->SplitAtomCardOfModelCard(chain->GetResidueChainId(), chain->GetEndingResidueSequenceNumber() + 1);
    }
    cout << "Putting TER card after non-amino acid residues: Done" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Putting TER card after non-amino acid residues: Done" );
    vector<string> glycam_residue_names = this->GetAllResidueNamesFromDatasetFiles(glycam_lib_files, prep_files);
    // Get all TER card positions and split
    vector<pair<char, int> > ter_card_positions = pdb_file->GetAllTerCardPositions(glycam_residue_names);
    for(vector<pair<char, int> >::iterator it1 = ter_card_positions.begin(); it1 != ter_card_positions.end(); it1++)
    {
        pair<char, int> ter_position = *it1;
        char chain_id = ter_position.first;
        int sequence_number = ter_position.second;
        pdb_file->SplitAtomCardOfModelCard(chain_id, sequence_number);
    }
    cout << "Putting TER card after residues with no tail or with at least 2 tails: Done" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Putting TER card after residues with no tail or with at least 2 tails: Done" );

    // Add Terminals
    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
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
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, chain->GetResidueChainId(),
                                                                      chain->GetEndingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
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
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, chain->GetResidueChainId(),
                                                                      chain->GetStartingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueBefore(pdb_atom_card);
                }
            }
            else
            {
                // Add n terminal residue at the beginning of the chain
                PossibleNChainTermination n_termination = chain->GetSelectedNTermination();
                string string_n_termination = chain->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, chain->GetResidueChainId(),
                                                                      chain->GetStartingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueBefore(pdb_atom_card);
                }
                // Add c terminal residue at the end of the chain
                PossibleCChainTermination c_termination = chain->GetSelectedCTermination();
                string string_c_termination = chain->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue_from_c_termination = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, amino_lib_files);
                if(lib_file_residue_from_c_termination != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms_from_c_termination = lib_file_residue_from_c_termination->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card_for_c_termination = new PdbAtomCard();
                    pdb_atom_card_for_c_termination->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map_for_c_termination;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector_for_c_termination = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it3 = lib_atoms_from_c_termination.begin(); it3 != lib_atoms_from_c_termination.end(); it3++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it3).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, chain->GetResidueChainId(),
                                                                      chain->GetEndingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map_for_c_termination[serial_number] = pdb_atom;
                        atom_vector_for_c_termination.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card_for_c_termination->SetAtoms(atom_map_for_c_termination);
                    pdb_atom_card_for_c_termination->SetOrderedAtoms(atom_vector_for_c_termination);
                    pdb_file->InsertResidueAfter(pdb_atom_card_for_c_termination);
                }
            }
        }
    }
    cout << "Add terminals: Done" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Add terminals: Done" );
}

void PdbPreprocessor::UpdateAminoAcidChainsWithTheGivenModelNumber(PdbFile *pdb_file, vector<string> amino_lib_files,
                                                                   vector<string> glycam_lib_files, vector<string> prep_files,
                                                                   PdbPreprocessorChainTerminationVector chain_terminations, int model_number)
{
    // Before non-amino-acid residue
    for(PdbPreprocessor::PdbPreprocessorChainTerminationVector::iterator it1 = chain_terminations.begin(); it1 != chain_terminations.end(); it1++)
    {
        PdbPreprocessorChainTermination* chain = (*it1);
        pdb_file->SplitAtomCardOfModelCard(chain->GetResidueChainId(), chain->GetEndingResidueSequenceNumber() + 1);
    }

    cout << "Putting TER card after non-amino acid residues: Done" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Putting TER card after non-amino acid residues: Done" );
    vector<string> glycam_residue_names = this->GetAllResidueNamesFromDatasetFiles(glycam_lib_files, prep_files);
    // Get all TER card positions and split
    vector<pair<char, int> > ter_card_positions = pdb_file->GetAllTerCardPositions(glycam_residue_names);

    for(vector<pair<char, int> >::iterator it1 = ter_card_positions.begin(); it1 != ter_card_positions.end(); it1++)
    {
        pair<char, int> ter_position = *it1;
        char chain_id = ter_position.first;
        int sequence_number = ter_position.second;
        pdb_file->SplitAtomCardOfModelCard(chain_id, sequence_number);
    }
    cout << "Putting TER card after residues with no tail or with at least 2 tails: Done" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Putting TER card after residues with no tail or with at least 2 tails: Done" );

    // Add Terminals
    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
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
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, chain->GetResidueChainId(),
                                                                      chain->GetEndingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueAfterWithTheGivenModelNumber(pdb_atom_card, model_number);
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
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, chain->GetResidueChainId(),
                                                                      chain->GetStartingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueBeforeWithTheGivenModelNumber(pdb_atom_card, model_number);
                }
            }
            else
            {
                // Add n terminal residue at the beginning of the chain
                PossibleNChainTermination n_termination = chain->GetSelectedNTermination();
                string string_n_termination = chain->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, chain->GetResidueChainId(),
                                                                      chain->GetStartingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueBeforeWithTheGivenModelNumber(pdb_atom_card, model_number);
                }
                // Add c terminal residue at the end of the chain
                PossibleCChainTermination c_termination = chain->GetSelectedCTermination();
                string string_c_termination = chain->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue_from_c_termination = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, amino_lib_files);
                if(lib_file_residue_from_c_termination != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms_from_c_termination = lib_file_residue_from_c_termination->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card_for_c_termination = new PdbAtomCard();
                    pdb_atom_card_for_c_termination->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map_for_c_termination;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector_for_c_termination = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it3 = lib_atoms_from_c_termination.begin(); it3 != lib_atoms_from_c_termination.end(); it3++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it3).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, chain->GetResidueChainId(),
                                                                      chain->GetEndingResidueSequenceNumber(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map_for_c_termination[serial_number] = pdb_atom;
                        atom_vector_for_c_termination.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card_for_c_termination->SetAtoms(atom_map_for_c_termination);
                    pdb_atom_card_for_c_termination->SetOrderedAtoms(atom_vector_for_c_termination);
                    pdb_file->InsertResidueAfterWithTheGivenModelNumber(pdb_atom_card_for_c_termination, model_number);
                }
            }
        }
    }
    cout << "Add terminals: Done" << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Add terminals: Done" );
}

bool PdbPreprocessor::ExtractGapsInAminoAcidChains(string pdb_file_path)
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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

bool PdbPreprocessor::ExtractGapsInAminoAcidChains(string pdb_file_path, vector<string> amino_lib_files)
{
    try
    {
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();

        PdbPreprocessorChainIdResidueMap all_chain_map_residue;
        PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
        PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
        {
            PdbFileSpace::PdbResidue* residue = *it;
            char chain_id = residue->GetResidueChainId();

            stringstream ss;
            ss << chain_id;
            all_chain_map_residue[ss.str()].push_back(residue);
        }

        // Slower version
//        vector<string> amino_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
        // Advanced version
        ResidueNameMap amino_lib_residue_names_map = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);
        amino_lib_residue_names_map["HIS"] = "HIS";
        PdbPreprocessorChainIdResidueMap all_chain_map_amino_acid_residue = PdbPreprocessorChainIdResidueMap();

        PdbFileSpace::PdbFile::PdbResidueVector residues;
        for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
        {
            int internal_amino_acid_chain_counter = 0;
            PdbPreprocessorChainIdResidueMap chain_map_residue = PdbPreprocessorChainIdResidueMap();
            string chain_id = (*it).first;
            residues = (*it).second;

            for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
            {
                PdbResidue* residue = (*it1);
                if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
                {
                    stringstream ss;
                    ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    all_chain_map_amino_acid_residue[chain_id].push_back(residue);
                    stringstream sss(chain_id);
                    char c_id;
                    sss >> c_id;
                    chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                    chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
                }
                else
                {
                    stringstream ss;
                    ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                    chain_map_residue[ss.str()].push_back(residue);
                    internal_amino_acid_chain_counter++;
                }
            }
            if(chain_map_residue.size() > 2)
            {
                cout << "There is an undefined protein in the middle of the chain" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is an undefined protein in the middle of the chain" );
                cout << "Pdb file is not processible at this time" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR, "Pdb file is not processible at this time" );


                return false;
            }
        }

        PdbFile::PdbResidueAtomsMap residue_atom_map;
        residue_atom_map = pdb_file->GetAllAtomsOfResidues();
        for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_amino_acid_residue.begin(); it != all_chain_map_amino_acid_residue.end(); it++)
        {
            string chain_id = (*it).first;
            PdbFileSpace::PdbFile::PdbResidueVector residues = (*it).second;
            char c_id;
            stringstream ss(chain_id);
            ss >> c_id;
            vector<int> sequence_numbers = chain_map_sequence_number[c_id];
            vector<char> insertion_codes = chain_map_insertion_code[c_id];
            vector<int>::iterator starting_sequence_number_iterator = min_element(sequence_numbers.begin(), sequence_numbers.end());
            vector<int>::iterator ending_sequence_number_iterator = max_element(sequence_numbers.begin(), sequence_numbers.end());

            for(PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
            {
                int dist = distance(residues.begin(), it1);
                if(dist != residues.size() - 1)
                {
                    PdbResidue* residue = (*it1);
                    PdbResidue* next_residue = *(it1 + 1);
                    int i = distance(residues.begin(), it1);
                    int j = distance(residues.begin(), it1 + 1);
                    PdbAtom* c_atom_of_residue = pdb_file->GetAtomOfResidueByName(residue, "C", residue_atom_map);
                    PdbAtom* n_atom_of_next_residue = pdb_file->GetAtomOfResidueByName(next_residue, "N", residue_atom_map);
                    double distance = INFINITY;
                    if(c_atom_of_residue != NULL && n_atom_of_next_residue != NULL)
                    {
                        GeometryTopology::Coordinate c_atom_coordinate = c_atom_of_residue->GetAtomOrthogonalCoordinate();
                        GeometryTopology::Coordinate n_atom_coordinate = n_atom_of_next_residue->GetAtomOrthogonalCoordinate();
                        distance = c_atom_coordinate.Distance(n_atom_coordinate);
                    }
                    if(distance > dCutOff + 1.0)
                    {
                        PdbPreprocessorMissingResidue* missing_residues = new PdbPreprocessorMissingResidue(c_id, *starting_sequence_number_iterator,
                                                                                                            *ending_sequence_number_iterator, sequence_numbers.at(i),
                                                                                                            sequence_numbers.at(j), insertion_codes.at(i), insertion_codes.at(j));
                        missing_residues_.push_back(missing_residues);
                    }
                }
            }
        }

        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

bool PdbPreprocessor::ExtractGapsInAminoAcidChains(PdbFile* pdb_file)
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
    return true;
}

bool PdbPreprocessor::ExtractGapsInAminoAcidChains(PdbFile *pdb_file, vector<string> amino_lib_files)
{
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResiduesFromAtomCard();

    PdbPreprocessorChainIdResidueMap all_chain_map_residue;
    PdbPreprocessor::PdbPreprocessorChainIdSequenceNumbersMap chain_map_sequence_number;
    PdbPreprocessor::PdbPreprocessorChainIdInsertionCodeMap chain_map_insertion_code;

    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbFileSpace::PdbResidue* residue = *it;
        char chain_id = residue->GetResidueChainId();

        stringstream ss;
        ss << chain_id;
        all_chain_map_residue[ss.str()].push_back(residue);
    }

    // Slower version
//        vector<string> amino_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
    // Advanced version
    ResidueNameMap amino_lib_residue_names_map = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);
    amino_lib_residue_names_map["HIS"] = "HIS";
    PdbPreprocessorChainIdResidueMap all_chain_map_amino_acid_residue = PdbPreprocessorChainIdResidueMap();

    PdbFileSpace::PdbFile::PdbResidueVector residues;
    PdbPreprocessorChainIdResidueMap chain_map_residue;
    for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_residue.begin(); it != all_chain_map_residue.end(); it++)
    {
        int internal_amino_acid_chain_counter = 0;
        chain_map_residue = PdbPreprocessorChainIdResidueMap();
        string chain_id = (*it).first;
        residues = (*it).second;        

        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            PdbResidue* residue = (*it1);
            if(amino_lib_residue_names_map.find(residue->GetResidueName()) != amino_lib_residue_names_map.end())
            {
                stringstream ss;
                ss << "A_" << chain_id << "_" << internal_amino_acid_chain_counter;
                chain_map_residue[ss.str()].push_back(residue);
                all_chain_map_amino_acid_residue[chain_id].push_back(residue);
                stringstream sss(chain_id);
                char c_id;
                sss >> c_id;
                chain_map_sequence_number[c_id].push_back(residue->GetResidueSequenceNumber());
                chain_map_insertion_code[c_id].push_back(residue->GetResidueInsertionCode());
            }
            else
            {
                stringstream ss;
                ss << "NA_" << chain_id << "_" << internal_amino_acid_chain_counter;
                chain_map_residue[ss.str()].push_back(residue);
                internal_amino_acid_chain_counter++;
            }
        }
        if(chain_map_residue.size() > 2)
        {
            cout << "There is an undefined protein in the middle of the chain" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "There is an undefined protein in the middle of the chain" );
            cout << "Pdb file is not processible at this time" << endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Pdb file is not processible at this time" );

            return false;
        }
    }


    PdbFile::PdbResidueAtomsMap residue_atom_map;
    residue_atom_map = pdb_file->GetAllAtomsOfResidues();
    for(PdbPreprocessorChainIdResidueMap::iterator it = all_chain_map_amino_acid_residue.begin(); it != all_chain_map_amino_acid_residue.end(); it++)
    {
        string chain_id = (*it).first;
        PdbFileSpace::PdbFile::PdbResidueVector residues = (*it).second;
        char c_id;
        stringstream ss(chain_id);
        ss >> c_id;
        vector<int> sequence_numbers = chain_map_sequence_number[c_id];
        vector<char> insertion_codes = chain_map_insertion_code[c_id];
        vector<int>::iterator starting_sequence_number_iterator = min_element(sequence_numbers.begin(), sequence_numbers.end());
        vector<int>::iterator ending_sequence_number_iterator = max_element(sequence_numbers.begin(), sequence_numbers.end());        

        for(PdbFile::PdbResidueVector::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
        {
            int dist = distance(residues.begin(), it1);
            if(dist != residues.size() - 1)
            {
                PdbResidue* residue = (*it1);
                PdbResidue* next_residue = *(it1 + 1);
                int i = distance(residues.begin(), it1);
                int j = distance(residues.begin(), it1 + 1);
                PdbAtom* c_atom_of_residue = pdb_file->GetAtomOfResidueByName(residue, "C", residue_atom_map);
                PdbAtom* n_atom_of_next_residue = pdb_file->GetAtomOfResidueByName(next_residue, "N", residue_atom_map);
                double distance = INFINITY;
                if(c_atom_of_residue != NULL && n_atom_of_next_residue != NULL)
                {
                    GeometryTopology::Coordinate c_atom_coordinate = c_atom_of_residue->GetAtomOrthogonalCoordinate();
                    GeometryTopology::Coordinate n_atom_coordinate = n_atom_of_next_residue->GetAtomOrthogonalCoordinate();
                    distance = c_atom_coordinate.Distance(n_atom_coordinate);
                }
                if(distance > dCutOff + 1.0)
                {
                    PdbPreprocessorMissingResidue* missing_residues = new PdbPreprocessorMissingResidue(c_id, *starting_sequence_number_iterator,
                                                                                                        *ending_sequence_number_iterator, sequence_numbers.at(i),
                                                                                                        sequence_numbers.at(j), insertion_codes.at(i), insertion_codes.at(j));
                    missing_residues_.push_back(missing_residues);
                }
            }
        }
    }

    return true;
}

void PdbPreprocessor::UpdateGapsInAminoAcidChains(PdbFile* pdb_file, vector<string> amino_lib_files, PdbPreprocessorMissingResidueVector gaps)
{
    for(PdbPreprocessor::PdbPreprocessorMissingResidueVector::iterator it1 = gaps.begin(); it1 != gaps.end(); it1++)
    {
        PdbPreprocessorMissingResidue* gap = (*it1);
        pdb_file->SplitAtomCardOfModelCard(gap->GetResidueChainId(), gap->GetResidueAfterGap());
    }
    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;    
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
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueBeforeGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
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
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueAfterGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
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
                // Add n terminal residue at the beginning of the chain
                PossibleNChainTermination n_termination = gap->GetSelectedNTermination();
                string string_n_termination = gap->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueAfterGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueBefore(pdb_atom_card);
                }
                // Add c terminal residue at the end of the chain
                PossibleCChainTermination c_termination = gap->GetSelectedCTermination();
                string string_c_termination = gap->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue_from_c_termination = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, amino_lib_files);
                if(lib_file_residue_from_c_termination != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms_from_c_termination = lib_file_residue_from_c_termination->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card_for_c_termination = new PdbAtomCard();
                    pdb_atom_card_for_c_termination->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map_for_c_termination;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector_for_c_termination = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it3 = lib_atoms_from_c_termination.begin(); it3 != lib_atoms_from_c_termination.end(); it3++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it3).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueBeforeGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map_for_c_termination[serial_number] = pdb_atom;
                        atom_vector_for_c_termination.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card_for_c_termination->SetAtoms(atom_map_for_c_termination);
                    pdb_atom_card_for_c_termination->SetOrderedAtoms(atom_vector_for_c_termination);
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
            }
        }
    }
}
void PdbPreprocessor::UpdateGapsInAminoAcidChainsWithTheGivenModelNumber(PdbFile *pdb_file, vector<string> amino_lib_files, PdbPreprocessorMissingResidueVector gaps, int model_number)
{
    for(PdbPreprocessor::PdbPreprocessorMissingResidueVector::iterator it1 = gaps.begin(); it1 != gaps.end(); it1++)
    {
        PdbPreprocessorMissingResidue* gap = (*it1);
        pdb_file->SplitAtomCardOfModelCardWithTheGivenModelNumber(gap->GetResidueChainId(), gap->GetResidueAfterGap(), model_number);
    }
    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
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
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueBeforeGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueAfterWithTheGivenModelNumber(pdb_atom_card, model_number);
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
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueAfterGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueBeforeWithTheGivenModelNumber(pdb_atom_card, model_number);
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
                PossibleNChainTermination n_termination = gap->GetSelectedNTermination();
                string string_n_termination = gap->GetStringFormatOfNTermination(n_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue = GetLibraryResidueByNameFromMultipleLibraryFiles(string_n_termination, amino_lib_files);
                if(lib_file_residue != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms = lib_file_residue->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card = new PdbAtomCard();
                    pdb_atom_card->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = lib_atoms.begin(); it2 != lib_atoms.end(); it2++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it2).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_n_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueAfterGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map[serial_number] = pdb_atom;
                        atom_vector.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card->SetAtoms(atom_map);
                    pdb_atom_card->SetOrderedAtoms(atom_vector);
                    pdb_file->InsertResidueBeforeWithTheGivenModelNumber(pdb_atom_card, model_number);
                }
                // Add c terminal residue at the end of the chain
                PossibleCChainTermination c_termination = gap->GetSelectedCTermination();
                string string_c_termination = gap->GetStringFormatOfCTermination(c_termination);
                LibraryFileSpace::LibraryFileResidue* lib_file_residue_from_c_termination = GetLibraryResidueByNameFromMultipleLibraryFiles(string_c_termination, amino_lib_files);
                if(lib_file_residue_from_c_termination != NULL)
                {
                    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms_from_c_termination = lib_file_residue_from_c_termination->GetAtoms();
                    PdbFileSpace::PdbAtomCard* pdb_atom_card_for_c_termination = new PdbAtomCard();
                    pdb_atom_card_for_c_termination->SetRecordName("ATOM");
                    PdbFileSpace::PdbAtomCard::PdbAtomMap atom_map_for_c_termination;
                    PdbFileSpace::PdbAtomCard::PdbAtomOrderVector atom_vector_for_c_termination = PdbFileSpace::PdbAtomCard::PdbAtomOrderVector();
                    int serial_number = 0;
                    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it3 = lib_atoms_from_c_termination.begin(); it3 != lib_atoms_from_c_termination.end(); it3++)
                    {
                        LibraryFileSpace::LibraryFileAtom* lib_file_atom = (*it3).second;
                        PdbFileSpace::PdbAtom* pdb_atom = new PdbAtom(serial_number, lib_file_atom->GetName(), ' ', string_c_termination, gap->GetResidueChainId(),
                                                                      gap->GetResidueBeforeGap(), ' ', lib_file_atom->GetCoordinate(), dNotSet, dNotSet, " ", "");
                        atom_map_for_c_termination[serial_number] = pdb_atom;
                        atom_vector_for_c_termination.push_back(pdb_atom);
                        serial_number++;
                    }
                    pdb_atom_card_for_c_termination->SetAtoms(atom_map_for_c_termination);
                    pdb_atom_card_for_c_termination->SetOrderedAtoms(atom_vector_for_c_termination);
                    pdb_file->InsertResidueAfterWithTheGivenModelNumber(pdb_atom_card_for_c_termination, model_number);
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
bool PdbPreprocessor::ExtractAlternateResidue(string pdb_file_path)
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
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractAlternateResidue(PdbFile* pdb_file)
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
    return true;
}
void PdbPreprocessor::RemoveUnselectedAlternateResidues(PdbFile *pdb_file, PdbPreprocessorAlternateResidueMap alternate_residue_map)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
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
//                pdb_file->DeleteResidue(pdb_residue);
                to_be_deleted_residues_.push_back(pdb_residue);
            }
        }
    }
    DeleteAllToBeDeletedEntities(pdb_file);
}
void PdbPreprocessor::RemoveUnselectedAlternateResiduesWithTheGivenModelNumber(PdbFile *pdb_file, PdbPreprocessorAlternateResidueMap alternate_residue_map, int model_number)
{
    to_be_deleted_atoms_.clear();
    to_be_deleted_residues_.clear();
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
//                pdb_file->DeleteResidueWithTheGivenModelNumber(pdb_residue, model_number);
                to_be_deleted_residues_.push_back(pdb_residue);
            }
        }
    }
    DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(pdb_file);
}

void PdbPreprocessor::DeleteAllToBeDeletedEntities(PdbFile *pdb_file)
{
    if(to_be_deleted_atoms_.size() != 0)
        pdb_file->DeleteAtoms(to_be_deleted_atoms_);

    if(to_be_deleted_residues_.size() != 0)
        pdb_file->DeleteResidues(to_be_deleted_residues_);
}
void PdbPreprocessor::DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(PdbFile *pdb_file, int model_number)
{
    if(to_be_deleted_atoms_.size() != 0)
        pdb_file->DeleteAtomsWithTheGivenModelNumber(to_be_deleted_atoms_, model_number);

    if(to_be_deleted_residues_.size() != 0)
        pdb_file->DeleteResiduesWithTheGivenModelNumber(to_be_deleted_residues_, model_number);
}
bool PdbPreprocessor::ExtractResidueInfo(string pdb_file_path, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    try
    {
        vector<string> lib_files = vector<string>();
        for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();
        LibraryFileSpace::LibraryFile::ResidueMap lib_residues = GetAllResiduesFromMultipleLibFilesMap(lib_files);
        PrepFileSpace::PrepFile::ResidueMap prep_residues = GetAllResiduesFromMultiplePrepFilesMap(prep_files);

        LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
        {
            double residue_charge = 0.0;
            PdbResidue* residue = *it;
            if(lib_residues.find(residue->GetResidueName()) != lib_residues.end())
            {
                LibraryFileSpace::LibraryFileResidue* lib_residue = lib_residues[residue->GetResidueName()];
    //            PdbFile::PdbAtomVector atoms = pdb_file->GetAllAtomsOfResidue(residue);

    //            for(PdbFile::PdbAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
    //            {
    //                PdbAtom* atom = (*it1);
    //                LibraryFileSpace::LibraryFileAtom* lib_atom = lib_residue->GetLibraryAtomByAtomName(atom->GetAtomName());
    //                if(lib_atom != NULL)
    //                {
    //                    if(lib_atom->GetCharge() != dNotSet)
    //                        residue_charge += lib_atom->GetCharge();
    //                }
    //            }
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
                PrepFileSpace::PrepFileResidue::PrepFileAtomVector prep_atoms = prep_residue->GetAtoms();
                for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
                {
                    PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
                    residue_charge += prep_atom->GetCharge();
                }
            }
            stringstream ss;
            ss << residue->GetResidueName() << "_" << residue->GetResidueChainId() << "_" << residue->GetResidueSequenceNumber() << "_" << residue->GetResidueInsertionCode()
               << "_" << residue->GetResidueAlternateLocation();
            PdbPreprocessorResidueInfo* residue_info = new PdbPreprocessorResidueInfo(residue->GetResidueName(), residue->GetResidueChainId(), residue->GetResidueSequenceNumber(),
                                                                                      residue->GetResidueInsertionCode(), residue->GetResidueAlternateLocation(), residue_charge);
            residue_info_map_[ss.str()] = residue_info;
        }
        return true;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
bool PdbPreprocessor::ExtractResidueInfo(PdbFile *pdb_file, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    vector<string> lib_files = vector<string>();
    for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    LibraryFileSpace::LibraryFile::ResidueMap lib_residues = GetAllResiduesFromMultipleLibFilesMap(lib_files);
    PrepFileSpace::PrepFile::ResidueMap prep_residues = GetAllResiduesFromMultiplePrepFilesMap(prep_files);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();

    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
    PrepFileSpace::PrepFileResidue::PrepFileAtomVector prep_atoms;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        double residue_charge = 0.0;
        PdbResidue* residue = *it;
        if(lib_residues.find(residue->GetResidueName()) != lib_residues.end())
        {
            LibraryFileSpace::LibraryFileResidue* lib_residue = lib_residues[residue->GetResidueName()];
//            PdbFile::PdbAtomVector atoms = pdb_file->GetAllAtomsOfResidue(residue);

//            for(PdbFile::PdbAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
//            {
//                PdbAtom* atom = (*it1);
//                LibraryFileSpace::LibraryFileAtom* lib_atom = lib_residue->GetLibraryAtomByAtomName(atom->GetAtomName());
//                if(lib_atom != NULL)
//                {
//                    if(lib_atom->GetCharge() != dNotSet)
//                        residue_charge += lib_atom->GetCharge();
//                }
//            }
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
            for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
            {
                PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
                residue_charge += prep_atom->GetCharge();
            }
        }
        stringstream ss;
        ss << residue->GetResidueName() << "_" << residue->GetResidueChainId() << "_" << residue->GetResidueSequenceNumber() << "_" << residue->GetResidueInsertionCode()
           << "_" << residue->GetResidueAlternateLocation();
        PdbPreprocessorResidueInfo* residue_info = new PdbPreprocessorResidueInfo(residue->GetResidueName(), residue->GetResidueChainId(), residue->GetResidueSequenceNumber(),
                                                                                  residue->GetResidueInsertionCode(), residue->GetResidueAlternateLocation(), residue_charge);
        residue_info_map_[ss.str()] = residue_info;
    }
    return true;
}
double PdbPreprocessor::CalculateModelCharge(string pdb_file_path, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    try
    {
        vector<string> lib_files = vector<string>();
        for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
        {
            lib_files.push_back(*it);
        }
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        double model_charge = 0.0;
        LibraryFileSpace::LibraryFile::ResidueMap lib_residues = GetAllResiduesFromMultipleLibFilesMap(lib_files);
        PrepFileSpace::PrepFile::ResidueMap prep_residues = GetAllResiduesFromMultiplePrepFilesMap(prep_files);
        PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();

        LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
        PrepFileSpace::PrepFileResidue::PrepFileAtomVector prep_atoms;
        for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
        {
            PdbResidue* residue = *it;
            if(lib_residues.find(residue->GetResidueName()) != lib_residues.end())
            {
                LibraryFileSpace::LibraryFileResidue* lib_residue = lib_residues[residue->GetResidueName()];
//                PdbFile::PdbAtomVector atoms = pdb_file->GetAllAtomsOfResidue(residue);

//                for(PdbFile::PdbAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
//                {
//                    PdbAtom* atom = (*it1);
//                    LibraryFileSpace::LibraryFileAtom* lib_atom = lib_residue->GetLibraryAtomByAtomName(atom->GetAtomName());
//                    if(lib_atom != NULL)
//                    {
//                        if(lib_atom->GetCharge() != dNotSet)
//                            model_charge += lib_atom->GetCharge();
//                    }
//                }

                lib_atoms = lib_residue->GetAtoms();
                for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it1 = lib_atoms.begin(); it1 != lib_atoms.end(); it1++)
                {
                    LibraryFileSpace::LibraryFileAtom* lib_atom = (*it1).second;
                    model_charge += lib_atom->GetCharge();
                }
            }
            else if (prep_residues.find(residue->GetResidueName()) != prep_residues.end())
            {
                PrepFileSpace::PrepFileResidue* prep_residue = prep_residues[residue->GetResidueName()];
                prep_atoms = prep_residue->GetAtoms();
                for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
                {
                    PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
                    model_charge += prep_atom->GetCharge();
                }
            }
        }
        return model_charge;
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}
double PdbPreprocessor::CalculateModelCharge(PdbFile* pdb_file, vector<string> amino_lib_files, vector<string> glycam_lib_files, vector<string> other_lib_files, vector<string> prep_files)
{
    vector<string> lib_files = vector<string>();
    for(vector<string>::iterator it = amino_lib_files.begin(); it != amino_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = glycam_lib_files.begin(); it != glycam_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    for(vector<string>::iterator it = other_lib_files.begin(); it != other_lib_files.end(); it++)
    {
        lib_files.push_back(*it);
    }
    double model_charge = 0.0;
    LibraryFileSpace::LibraryFile::ResidueMap lib_residues = GetAllResiduesFromMultipleLibFilesMap(lib_files);
    PrepFileSpace::PrepFile::ResidueMap prep_residues = GetAllResiduesFromMultiplePrepFilesMap(prep_files);
    PdbFileSpace::PdbFile::PdbResidueVector pdb_residues = pdb_file->GetAllResidues();

    LibraryFileSpace::LibraryFileResidue::AtomMap lib_atoms;
    PrepFileSpace::PrepFileResidue::PrepFileAtomVector prep_atoms;
    for(PdbFileSpace::PdbFile::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
    {
        PdbResidue* residue = *it;
        if(lib_residues.find(residue->GetResidueName()) != lib_residues.end())
        {
            LibraryFileSpace::LibraryFileResidue* lib_residue = lib_residues[residue->GetResidueName()];
//            PdbFile::PdbAtomVector atoms = pdb_file->GetAllAtomsOfResidue(residue);

//            for(PdbFile::PdbAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
//            {
//                PdbAtom* atom = (*it1);
//                LibraryFileSpace::LibraryFileAtom* lib_atom = lib_residue->GetLibraryAtomByAtomName(atom->GetAtomName());
//                if(lib_atom != NULL)
//                {
//                    if(lib_atom->GetCharge() != dNotSet)
//                        model_charge += lib_atom->GetCharge();
//                }
//            }

            lib_atoms = lib_residue->GetAtoms();
            for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it1 = lib_atoms.begin(); it1 != lib_atoms.end(); it1++)
            {
                LibraryFileSpace::LibraryFileAtom* lib_atom = (*it1).second;
                model_charge += lib_atom->GetCharge();
            }
        }
        else if (prep_residues.find(residue->GetResidueName()) != prep_residues.end())
        {
            PrepFileSpace::PrepFileResidue* prep_residue = prep_residues[residue->GetResidueName()];
            prep_atoms = prep_residue->GetAtoms();
            for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it1 = prep_atoms.begin(); it1 != prep_atoms.end(); it1++)
            {
                PrepFileSpace::PrepFileAtom* prep_atom = (*it1);
                model_charge += prep_atom->GetCharge();
            }
        }
    }
    return model_charge;
}

void PdbPreprocessor::Preprocess(PdbFile* pdb_file, vector<string> amino_lib_files_path, vector<string> glycam_lib_files_path, vector<string> other_lib_files_path, vector<string> prep_files_path)
{
    try
    {
        time_t t = time(0);
        string time_str = std::asctime(std::localtime(&t));
        cout << time_str.substr(0, time_str.size() - 1) << " Start preprocessing ..." << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, " Start preprocessing ..." );
        bool his_ext = ExtractHISResidues(pdb_file);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream his;
        if(his_ext)
            his << time_str.substr(0, time_str.size() - 1) << " HIS residues extraction: done";
        else
            his << time_str.substr(0, time_str.size() - 1) << " HIS residues extraction: failed";
        cout << his.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, his.str());
        bool cys_ext = ExtractCYSResidues(pdb_file);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream cys;
        if(cys_ext)
            cys << time_str.substr(0, time_str.size() - 1) << " CYS residues extraction: done";
        else
            cys << time_str.substr(0, time_str.size() - 1) << " CYS residues extraction: failed";
        cout << cys.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, cys.str() );
        bool alt_res_ext = ExtractAlternateResidue(pdb_file);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream alt;
        if(alt_res_ext)
            alt << time_str.substr(0, time_str.size() - 1) << " Alternate residues extraction: done";
        else
            alt << time_str.substr(0, time_str.size() - 1) << " Alternate residues extraction: failed";
        cout << alt.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, alt.str() );
        bool unrec_res_ext = ExtractUnrecognizedResidues(pdb_file, amino_lib_files_path, glycam_lib_files_path, other_lib_files_path, prep_files_path);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream unrecognized;
        if(unrec_res_ext)
            unrecognized << time_str.substr(0, time_str.size() - 1) << " Unrecognized residues extraction: done";
        else
            unrecognized << time_str.substr(0, time_str.size() - 1) << " Unrecognized residues extraction: failed";
        cout << unrecognized.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, unrecognized.str() );
        bool unknown_heavy_atom_ext = ExtractUnknownHeavyAtoms(pdb_file, amino_lib_files_path, glycam_lib_files_path, other_lib_files_path, prep_files_path);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream heavy;
        if(unknown_heavy_atom_ext)
            heavy << time_str.substr(0, time_str.size() - 1) << " Unknown heavy atoms extraction: done" ;
        else
            heavy << time_str.substr(0, time_str.size() - 1) << " Unknown heavy atoms extraction: failed" ;
        cout << heavy.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, heavy.str() );
        bool removed_hydro_ext = ExtractRemovedHydrogens(pdb_file, amino_lib_files_path, glycam_lib_files_path, other_lib_files_path, prep_files_path);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream hydrogen;
        if(removed_hydro_ext)
            hydrogen << time_str.substr(0, time_str.size() - 1) << " Removed hydrogens extraction: done" ;
        else
            hydrogen << time_str.substr(0, time_str.size() - 1) << " Removed hydrogens extraction: failed" ;
        cout << hydrogen.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, hydrogen.str() );
        bool amino_ext = ExtractAminoAcidChains(pdb_file, amino_lib_files_path);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream amino;
        if(amino_ext)
            amino << time_str.substr(0, time_str.size() - 1) << " Amino acid chains extraction: done" ;
        else
            amino << time_str.substr(0, time_str.size() - 1) << " Amino acid chains extraction: failed" ;
        cout << amino.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, amino.str() );
        bool gap_ext = ExtractGapsInAminoAcidChains(pdb_file, amino_lib_files_path);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream gaps;
        if(gap_ext)
            gaps << time_str.substr(0, time_str.size() - 1) << " Gaps in amino acid chains extraction: done";
        else
            gaps << time_str.substr(0, time_str.size() - 1) << " Gaps in amino acid chains extraction: failed";
        cout << gaps.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, gaps.str() );
        bool res_inf_ext = ExtractResidueInfo(pdb_file, amino_lib_files_path, glycam_lib_files_path, other_lib_files_path, prep_files_path);
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream info;
        if(res_inf_ext)
            info << time_str.substr(0, time_str.size() - 1) << " Residue info extraction: done" ;
        else
            info << time_str.substr(0, time_str.size() - 1) << " Residue info extraction: failed" ;
        cout << info.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, info.str() );
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream model_charge;
        model_charge << "Model charge is " << CalculateModelCharge(pdb_file, amino_lib_files_path, glycam_lib_files_path, other_lib_files_path, prep_files_path) ;
        cout << model_charge.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, model_charge.str() );
        stringstream model_done;
        model_done << time_str.substr(0, time_str.size() - 1) << " Model charge calculation: done" ;
        cout << model_done.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, model_done.str() );
        t = time(0);
        time_str = std::asctime(std::localtime(&t));
        stringstream pre;
        pre << time_str.substr(0, time_str.size() - 1) << " Preprocessing done";
        cout << pre.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, pre.str() );
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

void PdbPreprocessor::ApplyPreprocessing(PdbFile *pdb_file, vector<string> amino_lib_files_path, vector<string> glycam_lib_files_path, vector<string> prep_files_path)
{
    time_t t = time(0);
    string time_str = std::asctime(std::localtime(&t));
    stringstream changes;
    changes << time_str.substr(0, time_str.size() - 1) << " Start to apply changes ...";
    cout << changes.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, changes.str() );
    UpdateHISMapping(pdb_file,this->GetHistidineMappings());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream his_update;
    his_update << time_str.substr(0, time_str.size() - 1) << " HIS residues update: done" ;
    cout << his_update.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, his_update.str() );
    UpdateCYSResidues(pdb_file, this->GetDisulfideBonds());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream cys_update;
    cys_update << time_str.substr(0, time_str.size() - 1) << " CYS residues update: done" ;
    cout << cys_update.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, cys_update.str() );
    RemoveUnselectedAlternateResidues(pdb_file,this->GetAlternateResidueMap());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream alt_res;
    alt_res << time_str.substr(0, time_str.size() - 1) << " Unselected alternate residues removed: done" ;
    cout << alt_res.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, alt_res.str() );
    RemoveUnrecognizedResidues(pdb_file, this->GetUnrecognizedResidues());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream remove_res;
    remove_res << time_str.substr(0, time_str.size() - 1) << " Remove unrecognized residues: done" ;
    cout << remove_res.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, remove_res.str() );
    RemoveResiduesOfUnknownHeavyAtoms(pdb_file, this->GetUnrecognizedHeavyAtoms());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream remove_heavy;
    remove_heavy << time_str.substr(0, time_str.size() - 1) << " Unknown heavy atoms removed: done" ;
    cout << remove_heavy.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, remove_heavy.str() );
    RemoveRemovedHydrogens(pdb_file, this->GetReplacedHydrogens());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream remove_hydrogen;
    remove_hydrogen << time_str.substr(0, time_str.size() - 1) << " Removed hydrogens removed: done" ;
    cout << remove_hydrogen.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, remove_hydrogen.str() );
    UpdateAminoAcidChains(pdb_file, amino_lib_files_path, glycam_lib_files_path, prep_files_path, this->GetChainTerminations());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream amino_update;
    amino_update << time_str.substr(0, time_str.size() - 1) << " Amino acid chains update: done" ;
    cout << amino_update.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, amino_update.str() );
    UpdateGapsInAminoAcidChains(pdb_file, amino_lib_files_path, this->GetMissingResidues());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream gaps_update;
    gaps_update << time_str.substr(0, time_str.size() - 1) << " Gaps in amino acid chains update: done" ;
    cout << gaps_update.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, gaps_update.str() );
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream applied;
    applied << time_str.substr(0, time_str.size() - 1) << " Applying changes done" ;
    cout << applied.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, applied.str() );
}

void PdbPreprocessor::ApplyPreprocessingWithTheGivenModelNumber(PdbFile *pdb_file, vector<string> amino_lib_files_path,
                                                                vector<string> glycam_lib_files_path, vector<string> prep_files_path, int model_number)
{
    time_t t = time(0);
    string time_str = std::asctime(std::localtime(&t));
    stringstream changes;
    changes << time_str.substr(0, time_str.size() - 1) << " Start to apply changes ..." ;
    cout << changes.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, changes.str() );
    UpdateHISMappingWithTheGivenNumber(pdb_file,this->GetHistidineMappings(), model_number);
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream his_update;
    his_update << time_str.substr(0, time_str.size() - 1) << " HIS residues update: done" ;
    cout << his_update.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, his_update.str() );
    UpdateCYSResiduesWithTheGivenModelNumber(pdb_file, this->GetDisulfideBonds());
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream cys_update;
    cys_update << time_str.substr(0, time_str.size() - 1) << " CYS residues update: done" ;
    cout << cys_update.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, cys_update.str() );
    RemoveUnselectedAlternateResiduesWithTheGivenModelNumber(pdb_file,this->GetAlternateResidueMap(), model_number);
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream alt_res;
    alt_res << time_str.substr(0, time_str.size() - 1) << " Unselected alternate residues removed: done" ;
    cout << alt_res.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, alt_res.str() );
    RemoveUnrecognizedResiduesWithTheGivenModelNumber(pdb_file, this->GetUnrecognizedResidues(), model_number);
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream remove_res;
    remove_res << time_str.substr(0, time_str.size() - 1) << " Remove unrecognized residues: done" ;
    cout << remove_res.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, remove_res.str() );
    RemoveResiduesOfUnknownHeavyAtomsWithTheGivenModelNumber(pdb_file, this->GetUnrecognizedHeavyAtoms(), model_number);
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream remove_heavy;
    remove_heavy << time_str.substr(0, time_str.size() - 1) << " Unknown heavy atoms removed: done" ;
    cout << remove_heavy.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, remove_heavy.str() );
    RemoveRemovedHydrogensWithTheGivenModelNumber(pdb_file, this->GetReplacedHydrogens(), model_number);
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream remove_hydrogen;
    remove_hydrogen << time_str.substr(0, time_str.size() - 1) << " Removed hydrogens removed: done" ;
    cout << remove_hydrogen.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, remove_hydrogen.str() );
    UpdateAminoAcidChainsWithTheGivenModelNumber(pdb_file, amino_lib_files_path, glycam_lib_files_path, prep_files_path, this->GetChainTerminations(), model_number);
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream amino_update;
    amino_update << time_str.substr(0, time_str.size() - 1) << " Amino acid chains update: done" ;
    cout << amino_update.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, amino_update.str() );
    UpdateGapsInAminoAcidChainsWithTheGivenModelNumber(pdb_file, amino_lib_files_path, this->GetMissingResidues(), model_number);
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream gaps_update;
    gaps_update << time_str.substr(0, time_str.size() - 1) << " Gaps in amino acid chains update: done" ;
    cout << gaps_update.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, gaps_update.str() );
    t = time(0);
    time_str = std::asctime(std::localtime(&t));
    stringstream applied;
    applied << time_str.substr(0, time_str.size() - 1) << " Applying changes done" ;
    cout << applied.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, applied.str() );
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
    out << "================================ Residues Info ====================================" << endl;
    for(PdbPreprocessorResidueInfoMap::iterator it = residue_info_map_.begin(); it != residue_info_map_.end(); it++)
    {
        PdbPreprocessorResidueInfo* residue_info = (*it).second;
        residue_info->Print(out);
    }
}

