
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
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"

using namespace std;
using namespace PdbPreprocessorSpace;
using namespace PdbFileSpace;

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
PdbPreprocessor::PdbResidueVector PdbPreprocessor::GetUnrecognizedResidues(PdbPreprocessor::PdbResidueVector pdb_residues, vector<string> unrecognized_residue_names)
{
    PdbPreprocessor::PdbResidueVector unrecognized_residues;
    bool is_unrecognized = false;
    for(PdbPreprocessor::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
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
PdbPreprocessor::PdbResidueVector PdbPreprocessor::GetRecognizedResidues(PdbPreprocessor::PdbResidueVector pdb_residues, vector<string> recognized_residue_names)
{
    PdbPreprocessor::PdbResidueVector recognized_residues;
    bool is_recognized = false;
    for(PdbPreprocessor::PdbResidueVector::iterator it = pdb_residues.begin(); it != pdb_residues.end(); it++)
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

void PdbPreprocessor::ExtractUnrecognizedResidues(string pdb_file_path, vector<string> lib_files, vector<string> prep_files)
{
    vector<string> dataset_residue_names = GetAllResidueNamesFromDatasetFiles(lib_files, prep_files);
    PdbFile* pdb_file = new PdbFile(pdb_file_path);
    vector<string> pdb_residue_names = pdb_file->GetAllResidueNames();
    vector<string> unrecognized_residue_names = GetUnrecognizedResidueNames(pdb_residue_names, dataset_residue_names);
    vector<PdbResidue*> pdb_residues = pdb_file->GetAllResidues();
    vector<PdbResidue*> unrecognized_residues = GetUnrecognizedResidues(pdb_residues, unrecognized_residue_names);

    for(vector<PdbResidue*>::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
    {
        PdbResidue* pdb_residue = (*it);
        PdbPreprocessorUnrecognizedResidue* unrecognized_residue =
                new PdbPreprocessorUnrecognizedResidue(pdb_residue->GetResidueName(), pdb_residue->GetResidueChainId(), pdb_residue->GetResidueSequenceNumber());
        unrecognized_residues_.push_back(unrecognized_residue);
    }
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
    return all_lib_residue_names;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbPreprocessor::Print(ostream &out)
{
}

