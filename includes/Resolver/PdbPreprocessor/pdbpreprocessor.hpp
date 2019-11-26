#ifndef PDBPREPROCESSOR_HPP
#define PDBPREPROCESSOR_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../../InputSet//PdbFileSpace/pdbresidue.hpp"
#include "../../InputSet//PdbFileSpace/pdbfile.hpp"
#include "../../InputSet//PdbFileSpace/pdbatomcard.hpp"
#include "../../ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../common.hpp"

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorDisulfideBond;
    class PdbPreprocessorChainTermination;
    class PdbPreprocessorHistidineMapping;
    class PdbPreprocessorMissingResidue;
    class PdbPreprocessorUnrecognizedResidue;
    class PdbPreprocessorUnrecognizedHeavyAtom;
    class PdbPreprocessorReplacedHydrogen;
    class PdbPreprocessorAlternateResidue;
    class PdbPreprocessorResidueInfo;
    class PdbPreprocessor
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of disulfide bonds
              */
            typedef std::vector<PdbPreprocessorDisulfideBond*> PdbPreprocessorDisulfideBondVector;
            /*! \typedef
              * List of chain terminations
              */
            typedef std::vector<PdbPreprocessorChainTermination*> PdbPreprocessorChainTerminationVector;
            /*! \typedef
              * List of histidine residues
              */
            typedef std::vector<PdbPreprocessorHistidineMapping*> PdbPreprocessorHistidineMappingVector;
            /*! \typedef
              * List of missing residues (gaps)
              */
            typedef std::vector<PdbPreprocessorMissingResidue*> PdbPreprocessorMissingResidueVector;
            /*! \typedef
              * List of unrecognized residues
              */
            typedef std::vector<PdbPreprocessorUnrecognizedResidue*> PdbPreprocessorUnrecognizedResidueVector;
            /*! \typedef
              * List of Recognized residues
              */
            typedef std::vector<PdbPreprocessorUnrecognizedResidue*> PdbPreprocessorRecognizedResidueVector;
            /*! \typedef
              * List of unrecognized heavy atoms
              */
            typedef std::vector<PdbPreprocessorUnrecognizedHeavyAtom*> PdbPreprocessorUnrecognizedHeavyAtomVector;
            /*! \typedef
              * List of replaced/removed hydrogen atoms
              */
            typedef std::vector<PdbPreprocessorReplacedHydrogen*> PdbPreprocessorReplacedHydrogenVector;
            /*! \typedef
              * A mapping between residue chain id and its sequence number in that chain
              */
            typedef std::map<char, std::vector<int> > PdbPreprocessorChainIdSequenceNumbersMap;
            /*! \typedef
              * A mapping between residue chian id and its insertion code in that chain
              */
            typedef std::map<char, std::vector<char> > PdbPreprocessorChainIdInsertionCodeMap;
            /*! \typedef
              * A mapping between residue key with alternate location(s) and its corresponding alternate residue object(s)
              */
            typedef std::map<std::string, PdbPreprocessorAlternateResidue*> PdbPreprocessorAlternateResidueMap;
            /*! \typedef
              * List of atoms to be deleted
              */
            typedef std::vector<PdbFileSpace::PdbAtomCard*> PdbPreprocessorToBeDeletedAtomVector;
            /*! \typedef
              * List of residues to be deleted
              */
            typedef std::vector<PdbFileSpace::PdbResidue*> PdbPreprocessorToBeDeletedResidueVector;
            /*! \typedef
              * A mapping between residue chain id and the residues involving in the chain
              */
            typedef std::map<std::string, PdbFileSpace::PdbFile::PdbResidueVector > PdbPreprocessorChainIdResidueMap;
            /*! \typedef
              * A mapping between residue key and its corresponding residue info object
              */
            typedef std::map<std::string, PdbPreprocessorResidueInfo*> PdbPreprocessorResidueInfoMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessor();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the disulfide bonds
              * @return disulfide_bonds_ attribute of the current object of this class
              */
            PdbPreprocessorDisulfideBondVector GetDisulfideBonds();
            /*! \fn
              * An accessor function in order to access to the chain terminations
              * @return chain_terminations_ attribute of the current object of this class
              */
            PdbPreprocessorChainTerminationVector GetChainTerminations();
            /*! \fn
              * An accessor function in order to access to the histidine mappings
              * @return histidine_mappings_ attribute of the current object of this class
              */
            PdbPreprocessorHistidineMappingVector GetHistidineMappings();
            /*! \fn
              * An accessor function in order to access to the missing residues
              * @return missing_residues_ attribute of the current object of this class
              */
            PdbPreprocessorMissingResidueVector GetMissingResidues();
            /*! \fn
              * An accessor function in order to access to the unrecognized residues
              * @return unrecognized_residues_ attribute of the current object of this class
              */
            PdbPreprocessorUnrecognizedResidueVector GetUnrecognizedResidues();
            /*! \fn
              * An accessor function in order to access to the recognized residues
              * @return recognized_residues_ attribute of the current object of this class
              */
            PdbPreprocessorRecognizedResidueVector GetRecognizedResidues();
            /*! \fn
              * An accessor function in order to access to the unrecognized heavy atoms
              * @return unrecognized_heavy_atoms_ attribute of the current object of this class
              */
            PdbPreprocessorUnrecognizedHeavyAtomVector GetUnrecognizedHeavyAtoms();
            /*! \fn
              * An accessor function in order to access to the replaced hydrogens
              * @return replaced_hydrogens_ attribute of the current object of this class
              */
            PdbPreprocessorReplacedHydrogenVector GetReplacedHydrogens();
            /*! \fn
              * An accessor function in order to access to the alternate residue map
              * @return alternate_residue_map_ attribute of the current object of this class
              */
            PdbPreprocessorAlternateResidueMap GetAlternateResidueMap();
            /*! \fn
              * An accessor function in order to access to the list of to be deleted atoms
              * @return to_be_deleted_atoms_ attribute of the current object of this class
              */
            PdbPreprocessorToBeDeletedAtomVector GetToBeDeletedAtoms();
            /*! \fn
              * An accessor function in order to access to the list of to be deleted residues
              * @return to_be_deleted_residues_ attribute of the current object of this class
              */
            PdbPreprocessorToBeDeletedResidueVector GetToBeDeletedResidues();
            /*! \fn
              * An accessor function in order to access to the residue info map
              * @return reisude_info_map_ attribute of the current object of this class
              */
            PdbPreprocessorResidueInfoMap GetResidueInfoMap();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the disulfide bonds of the current object
              * Set the disulfide_bonds_ attribute of the current pdb preprocessor
              * @param disulfide_bonds The disulfide bonds attribute of the current object
              */
            void SetDisulfideBonds(PdbPreprocessorDisulfideBondVector disulfide_bonds);
            /*! \fn
              * A function in order to add the disulfide_bond to the current object
              * Set disulfide_bonds_ attribute of the current pdb preprocessor
              * @param disulfide_bond The heterogen atom attribute of the current object
              */
            void AddDisulfideBond(PdbPreprocessorDisulfideBond* disulfide_bond);
            /*! \fn
              * A mutator function in order to set the chain terminations of the current object
              * Set the chain_terminations_ attribute of the current pdb preprocessor
              * @param chain_terminations The chain terminations attribute of the current object
              */
            void SetChainTerminations(PdbPreprocessorChainTerminationVector chain_terminations);
            /*! \fn
              * A function in order to add the chain_termination to the current object
              * Set chain_terminations_ attribute of the current pdb preprocessor
              * @param chain_termination The chain termination attribute of the current object
              */
            void AddChainTermination(PdbPreprocessorChainTermination* chain_termination);
            /*! \fn
              * A mutator function in order to set the histidnie mappings of the current object
              * Set the histidnie_mappings_ attribute of the current pdb preprocessor
              * @param histidnie_mappings The histidnie mappings attribute of the current object
              */
            void SetHistidineMappings(PdbPreprocessorHistidineMappingVector histidine_mappings);
            /*! \fn
              * A function in order to add the histidnie mapping to the current object
              * Set histidnie_mappings_ attribute of the current pdb preprocessor
              * @param histidnie_mapping The histidnie mapping attribute of the current object
              */
            void AddHistidineMapping(PdbPreprocessorHistidineMapping* histidnine_mapping);
            /*! \fn
              * A mutator function in order to set the missing residues of the current object
              * Set the missing_residues_ attribute of the current pdb preprocessor
              * @param missing_residues The histidnie mappings attribute of the current object
              */
            void SetMissingResidues(PdbPreprocessorMissingResidueVector missing_residues);
            /*! \fn
              * A function in order to add the missing residue to the current object
              * Set missing_residues_ attribute of the current pdb preprocessor
              * @param missing_residue The histidnie mapping attribute of the current object
              */
            void AddMissingResidue(PdbPreprocessorMissingResidue* missing_residue);
            /*! \fn
              * A mutator function in order to set the unrecognized residues of the current object
              * Set the unrecognized_residues_ attribute of the current pdb preprocessor
              * @param unrecognized_residues The unrecognized residues attribute of the current object
              */
            void SetUnrecognizedResidues(PdbPreprocessorUnrecognizedResidueVector unrecognized_residues);
            /*! \fn
              * A function in order to add the unrecognized residue to the current object
              * Set unrecognized_residues_ attribute of the current pdb preprocessor
              * @param unrecognized_residue The unrecognized residue attribute of the current object
              */
            void AddUnrecognizedResidue(PdbPreprocessorUnrecognizedResidue* unrecognized_residue);
            /*! \fn
              * A mutator function in order to set the recognized residues of the current object
              * Set the recognized_residues_ attribute of the current pdb preprocessor
              * @param recognized_residues The recognized residues attribute of the current object
              */
            void SetRecognizedResidues(PdbPreprocessorRecognizedResidueVector recognized_residues);
            /*! \fn
              * A function in order to add the recognized residue to the current object
              * Set recognized_residues_ attribute of the current pdb preprocessor
              * @param recognized_residue The recognized residue attribute of the current object
              */
            void AddRecognizedResidue(PdbPreprocessorUnrecognizedResidue* recognized_residue);
            /*! \fn
              * A mutator function in order to set the unrecognized heavy atoms of the current object
              * Set the unrecognized_heavy_atoms_ attribute of the current pdb preprocessor
              * @param unrecognized_heavy_atoms The unrecognized heavy atoms attribute of the current object
              */
            void SetUnrecognizedHeavyAtoms(PdbPreprocessorUnrecognizedHeavyAtomVector unrecognized_heavy_atoms);
            /*! \fn
              * A function in order to add the unrecognized heavy atom to the current object
              * Set unrecognized_heavy_atoms_ attribute of the current pdb preprocessor
              * @param unrecognized_heavy_atom The unrecognized heavy atom attribute of the current object
              */
            void AddUnrecognizedHeavyAtom(PdbPreprocessorUnrecognizedHeavyAtom* unrecognized_heavy_atom);
            /*! \fn
              * A mutator function in order to set the replaced hydrogens of the current object
              * Set the replaced_hydrogens_ attribute of the current pdb preprocessor
              * @param replaced_hydrogens The replaced hydrogens attribute of the current object
              */
            void SetReplacedHydrogens(PdbPreprocessorReplacedHydrogenVector replaced_hydrogens);
            /*! \fn
              * A function in order to add the replaced hydrogen to the current object
              * Set replaced_hydrogens_ attribute of the current pdb preprocessor
              * @param replaced_hydrogen The replaced hydrogen atom attribute of the current object
              */
            void AddReplacedHydrogen(PdbPreprocessorReplacedHydrogen* replaced_hydrogen);
            /*! \fn
              * A mutator function in order to set the to be deleted atoms of the current object
              * Set the to_be_deleted_atoms attribute of the current pdb preprocessor
              * @param to_be_deleted_atoms The to be deleted atoms attribute of the current object
              */
            void SetToBeDeletedAtoms(PdbPreprocessorToBeDeletedAtomVector to_be_deleted_atoms);
            /*! \fn
              * A mutator function in order to set the to be deleted residues of the current object
              * Set the to_be_deleted_residues attribute of the current pdb preprocessor
              * @param to_be_deleted_residues The to be deleted residues attribute of the current object
              */
            void SetToBeDeletedResidues(PdbPreprocessorToBeDeletedResidueVector to_be_deleted_residues);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * A function in order to access to the unrecognized residue names of pdb
              * @param pdb_residue_names The list of residue names in the current pdb file
              * @param data_set_residue_names The list of residue names from library and prep files
              * @return unrecognized_residue_names
              */
            std::vector<std::string> GetUnrecognizedResidueNames(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, std::vector<std::string> dataset_residue_names);

            //**************************************************
            gmml::ResidueNameMap GetUnrecognizedResidueNamesMap(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, gmml::ResidueNameMap dataset_residue_names);

            //**************************************************

            /*! \fn
              * A function in order to access to the recognized residue names of pdb
              * @param pdb_residue_names The list of residue names in the current pdb file
              * @param data_set_residue_names The list of residue names from library and prep files
              * @return recognized_residue_names
              */
            std::vector<std::string> GetRecognizedResidueNames(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, std::vector<std::string> dataset_residue_names);

            //**************************************************
            gmml::ResidueNameMap GetRecognizedResidueNamesMap(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names, gmml::ResidueNameMap dataset_residue_names);

            //**************************************************

            /*! \fn
              * A function in order to access to the unrecognized residues of pdb
              * @param pdb_residues The list of residues in the current pdb file
              * @param unrecognized_residue_names The list of unrecognized residue names
              * @return unrecognized_residues
              */
            PdbFileSpace::PdbFile::PdbResidueVector GetUnrecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, std::vector<std::string> unrecognized_residue_names);

            //**************************************************
            PdbFileSpace::PdbFile::PdbResidueVector GetUnrecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, gmml::ResidueNameMap unrecognized_residue_names);

            //**************************************************

            /*! \fn
              * A function in order to access to the recognized residues of pdb
              * @param pdb_residues The list of residues in the current pdb file
              * @param recognized_residue_names The list of recognized residue names
              * @return recognized_residues
              */
            PdbFileSpace::PdbFile::PdbResidueVector GetRecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, std::vector<std::string> recognized_residue_names);

            //**************************************************
            PdbFileSpace::PdbFile::PdbResidueVector GetRecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, gmml::ResidueNameMap recognized_residue_names);

            //**************************************************

            /*! \fn
              * A function in order to access to the list of all residue names from lib files
              * @param lib_files The list of paths to library files
              * @return all_residue_names
              */
            std::vector<std::string> GetAllResidueNamesFromMultipleLibFiles(std::vector<std::string> lib_files);

            //**************************************************
            gmml::ResidueNameMap GetAllResidueNamesFromMultipleLibFilesMap(std::vector<std::string> lib_files);

            //**************************************************

            /*! \fn
              * A function in order to access to the list of all residue from lib files
              * @param lib_files The list of paths to library files
              * @return all_residues_
              */
            LibraryFileSpace::LibraryFile::ResidueMap GetAllResiduesFromMultipleLibFilesMap(std::vector<std::string> lib_files);
            /*! \fn
              * A function in order to access to the list of all residue from prep files
              * @param prep_files The list of paths to prep files
              * @return all_residues_
              */
            PrepFileSpace::PrepFile::ResidueMap GetAllResiduesFromMultiplePrepFilesMap(std::vector<std::string> prep_files);

            /*! \fn
              * A function in order to access to the list of all residue names from prep files
              * @param prep_files The list of paths to prep files
              * @return all_residue_names
              */
            std::vector<std::string> GetAllResidueNamesFromMultiplePrepFiles(std::vector<std::string> prep_files);

            //**************************************************
            gmml::ResidueNameMap GetAllResidueNamesFromMultiplePrepFilesMap(std::vector<std::string> prep_files);

            //**************************************************

            /*! \fn
              * A function in order to access to the list of all residue names from library and prep files
              * @param lib_files The list of paths to library files
              * @param prep_files The list of paths to prep files
              * @return all_residue_names
              */
            std::vector<std::string> GetAllResidueNamesFromDatasetFiles(std::vector<std::string> lib_files, std::vector<std::string> prep_files);

            //**************************************************
            gmml::ResidueNameMap GetAllResidueNamesFromDatasetFilesMap(std::vector<std::string> lib_files, std::vector<std::string> prep_files);

            //**************************************************
/** @}*/
            /** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A function in order to extract the unrecognized residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              * @return bool value
              */
            bool ExtractUnrecognizedResidues(std::string pdb_file_path, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to extract the unrecognized residues of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              * @return bool value
              */
            bool ExtractUnrecognizedResidues(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to remove the unrecognized residues of a pdb file
              * @param pdb_file The object of a pdb file
              * @param unknown_heavy_atoms The list of unrecognized residues
              */
            void RemoveUnrecognizedResidues(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorUnrecognizedResidueVector unrecognized_residues);
            /*! \fn
              * A function in order to remove the unrecognized residues of a pdb file
              * @param pdb_file The object of a pdb file
              * @param unknown_heavy_atoms The list of unrecognized residues
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void RemoveUnrecognizedResiduesWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorUnrecognizedResidueVector unrecognized_residues, int model_number = 1);
            /*! \fn
              * A function in order to extract the recognized residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              * @return bool value
              */
            bool ExtractRecognizedResidues(std::string pdb_file_path, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to extract the recognized residues of a pdb file
              * @param pdb_file The object of pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              * @return bool value
              */
            bool ExtractRecognizedResidues(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to access to the list of CYS residues
              * @param pdb_residues The list of pdb residues
              * @return all_cys_residues
              */
/** @}*/
            /** \addtogroup Molecular_Data_Structure
               * @{
               */
            PdbFileSpace::PdbFile::PdbResidueVector GetAllCYSResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues);
            /*! \fn
              * A function in order to access to the distance of a pair of CYS residues
              * @param first_residue The first residue of CYS pair
              * @param second_residue The second residue of CYS pair
              * @param pdb_file Pdb file object
              * @param residue_atom_map A map between a residue and its atoms
              * @param first_sulfur_atom_serial_number Serial number of the sulfur atom in the first residue
              * @param second_sulfur_atom_serial_number Serial number of the sulfur atom in the second residue
              * @return distance
              */
            double GetDistanceofCYS(PdbFileSpace::PdbResidue* first_residue, PdbFileSpace::PdbResidue* second_residue,
                                    PdbFileSpace::PdbFile* pdb_file, PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atom_map,
                                    int& first_sulfur_atom_serial_number, int& second_sulfur_atom_serial_number);
            /*! \fn
              * A function in order to extract the CYS residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @return bool value
              */
/** @}*/
            /** \addtogroup Manipulators
               * @{
               */

            bool ExtractCYSResidues(std::string pdb_file_path);
            /*! \fn
              * A function in order to extract the CYS residues of a pdb file
              * @param pdb_file The object of pdb file
              * @return bool value
              */
            bool ExtractCYSResidues(PdbFileSpace::PdbFile* pdb_file);
            /*! \fn
              * A function in order to remove the unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param disulfide_bonds The list of disulfide bonds
              */
            void UpdateCYSResidues(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorDisulfideBondVector disulfide_bonds);
            /*! \fn
             * A function in order to remove the unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param disulfide_bonds The list of disulfide bonds
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void UpdateCYSResiduesWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorDisulfideBondVector disulfide_bonds, int model_number = 1);
            /*! \fn
              * A function in order to access to the list of HIS residues
              * @param pdb_residues The list of pdb residues
              * @return all_his_residues
              */
/** @}*/
            /** \addtogroup Molecular_Data_Structure
               * @{
               */
            PdbFileSpace::PdbFile::PdbResidueVector GetAllHISResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues);
            /*! \fn
              * A function in order to extract the HIS residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @return bool value
              */
/** @}*/
            /** \addtogroup Manipulators
               * @{
               */
            bool ExtractHISResidues(std::string pdb_file_path);
            /*! \fn
              * A function in order to extract the HIS residues of a pdb file
              * @param pdb_file The object of pdb file
              * @return bool value
              */
            bool ExtractHISResidues(PdbFileSpace::PdbFile* pdb_file);
            /*! \fn
              * A function in order to update histidine mapping of a pdb file
              * @param pdb_file The object of a pdb file
              * @param histidine_mappings The list of histidine mappings
              */
            void UpdateHISMapping(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorHistidineMappingVector histidine_mappings);
            /*! \fn
              * A function in order to update histidine mapping of a pdb file
              * @param pdb_file The object of a pdb file
              * @param histidine_mappings The list of histidine mappings
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void UpdateHISMappingWithTheGivenNumber(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorHistidineMappingVector histidine_mappings, int model_number = 1);
            /*! \fn
              * A function in order to access to the list of unknown heavy atoms of a residue
              * @param pdb_atom_names_of_residue The list of atom names of a pdb residue
              * @param dataset_atom_names_of_residue The list of atom names of a residue from dataset
              * @return unknown_heavy_atom_names_of_residue
              */
/** @}*/
            /** \addtogroup Molecular_Data_Structure
               * @{
               */
            std::vector<std::string> GetUnknownHeavyAtomNamesOfResidue(std::vector<std::string> pdb_atom_names_of_residue, std::vector<std::string> dataset_atom_names_of_residue);

            /*! \fn
              * A function in order to access to the list of atom names of a residue from multiple library files
              * @param residue_name The name of the residue
              * @param lib_files The list of paths to the library files
              * @return all_atom_names_of_residue
              */
            std::vector<std::string> GetAllAtomNamesOfResidueFromMultipleLibFiles(std::string residue_name, std::vector<std::string> lib_files);

            //**************************************************
            gmml::ResidueNameAtomNamesMap GetAllAtomNamesOfResidueNamesFromMultipleLibFiles(std::vector<std::string> lib_files);
            //**************************************************

            /*! \fn
              * A function in order to access to the list of atom names of a residue from multiple prep files
              * @param residue_name The name of the residue
              * @param prep_files The list of paths to the prep files
              * @return all_atom_names_of_residue
              */
            std::vector<std::string> GetAllAtomNamesOfResidueFromMultiplePrepFiles(std::string residue_name, std::vector<std::string> prep_files);

            //**************************************************
            gmml::ResidueNameAtomNamesMap GetAllAtomNamesOfResidueNamesFromMultiplePrepFiles(std::vector<std::string> prep_files);
            //**************************************************

            /*! \fn
              * A function in order to access to the list of all atom names of a residue from library and prep files
              * @param lib_files The list of paths to library files
              * @param prep_files The list of paths to prep files
              * @return all_atom_names
              */
            std::vector<std::string> GetAllAtomNamesOfResidueFromDatasetFiles(std::string residue_name, std::vector<std::string> lib_files, std::vector<std::string> prep_files);

            //**************************************************
            gmml::ResidueNameAtomNamesMap GetAllAtomNamesOfResidueNamesFromDatasetFiles(std::vector<std::string> lib_files, std::vector<std::string> prep_files);
            //**************************************************

            /*! \fn
              * A function in order to access to the unknown heavy atoms of a residue
              * @param pdb_atoms The list of pbd atoms
              * @param dataset_atom_names_of_residue The list of atom names of a residue from dataset
              * @return unknown_heavy_atoms_of_residue
              */
            PdbFileSpace::PdbFile::PdbAtomCardVector GetUnknownHeavyAtomsOfResidue(PdbFileSpace::PdbFile::PdbAtomCardVector pdb_atoms, std::vector<std::string> dataset_atom_names_of_residue);

            /*! \fn
              * A function in order to extract the unknown heavy atoms of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              * @return bool value
              */
/** @}*/
            /** \addtogroup Manipulators
               * @{
               */
            bool ExtractUnknownHeavyAtoms(std::string pdb_file_path, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to extract the unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              * @return bool value
              */
            bool ExtractUnknownHeavyAtoms(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to remove the unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param unknown_heavy_atoms The list of unknown heavy atoms
              */
            void RemoveUnknownHeavyAtoms(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms);
            /*! \fn
              * A function in order to remove the unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param unknown_heavy_atoms The list of unknown heavy atoms
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void RemoveUnknownHeavyAtomsWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms, int model_number = 1);
            /*! \fn
              * A function in order to remove residues of unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param unknown_heavy_atoms The list of unknown heavy atoms
              */
            void RemoveResiduesOfUnknownHeavyAtoms(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms);
            /*! \fn
              * A function in order to remove residues of unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param unknown_heavy_atoms The list of unknown heavy atoms
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void RemoveResiduesOfUnknownHeavyAtomsWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms, int model_number = 1);
            /*! \fn
              * A function in order to access to the removed hydrogen names of a residue
              * @param pdb_atom_names_of_residue The list of atom names of a residue from pdb file
              * @param dataset_atom_names_of_residue The list of atom names of a residue from dataset
              * @return removed_hydrogen_names_of_residue
              */
/** @}*/
            /** \addtogroup Molecular_Data_Structure
               * @{
               */
            std::vector<std::string> GetRemovedHydrogenNamesOfResidue(std::vector<std::string> pdb_atom_names_of_residue, std::vector<std::string> dataset_atom_names_of_residue);

            /*! \fn
              * A function in order to access to the removed hydrogens of a residue
              * @param pdb_atoms The list of pdb atoms
              * @param dataset_atom_names_of_residue The list of atom names of a residue from dataset
              * @return removed_hydrogens_of_residue
              */
            PdbFileSpace::PdbFile::PdbAtomCardVector GetRemovedHydrogensOfResidue(PdbFileSpace::PdbFile::PdbAtomCardVector pdb_atoms, std::vector<std::string> dataset_atom_names_of_residue);

            /*! \fn
              * A function in order to extract the removed hydrogens of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              * @return bool value
              */
/** @}*/
            /** \addtogroup Manipulators
               * @{
               */
            bool ExtractRemovedHydrogens(std::string pdb_file_path, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to extract the removed hydrogens of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              * @return bool value
              */
            bool ExtractRemovedHydrogens(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to remove the removed hydrogens of a pdb file
              * @param pdb_file The object of a pdb file
              * @param replaced_hydrogens The list of replaced hydrogens
              */
            void RemoveRemovedHydrogens(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorReplacedHydrogenVector replaced_hydrogens);
            /*! \fn
              * A function in order to remove the removed hydrogens of a pdb file
              * @param pdb_file The object of a pdb file
              * @param replaced_hydrogens The list of replaced hydrogens
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void RemoveRemovedHydrogensWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorReplacedHydrogenVector replaced_hydrogens, int model_number = 1);
            /*! \fn
              * A function in order to extract the amino acid chains of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @return bool value
              */
            bool ExtractAminoAcidChains(std::string pdb_file_path);

            //**************************************************
            bool ExtractAminoAcidChains(std::string pdb_file_path, std::vector<std::string> amino_lib_files);

            //**************************************************
            /*! \fn
              * A function in order to extract the amino acid chains of a pdb file
              * @param pdb_file The object of pdb file
              * @return bool value
              */
            bool ExtractAminoAcidChains(PdbFileSpace::PdbFile* pdb_file);

            //**************************************************
            bool ExtractAminoAcidChains(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files);

            //**************************************************
            /*! \fn
              * A function in order to update the amino acid chains of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param gaps The list of chain terminations
              */
            void UpdateAminoAcidChains(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files,
                                       std::vector<std::string> prep_files, PdbPreprocessorChainTerminationVector chain_terminations);
            /*! \fn
              * A function in order to update the amino acid chains of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param gaps The list of chain terminations
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void UpdateAminoAcidChainsWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files,
                                                              std::vector<std::string> glycam_lib_files, std::vector<std::string> prep_files,
                                                              PdbPreprocessorChainTerminationVector chain_terminations, int model_number = 1);
            /*! \fn
              * A function in order to extract the gaps in amino acid chains of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @return bool value
              */
            bool ExtractGapsInAminoAcidChains(std::string pdb_file_path);

            //***************************************************
             bool ExtractGapsInAminoAcidChains(std::string pdb_file_path, std::vector<std::string> amino_lib_files);

            //************************************************
            /*! \fn
              * A function in order to extract the gaps in amino acid chains of a pdb file
              * @param pdb_file The object of pdb file
              * @return bool value
              */
            bool ExtractGapsInAminoAcidChains(PdbFileSpace::PdbFile* pdb_file);

            //***************************************************
             bool ExtractGapsInAminoAcidChains(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files);

            //************************************************
            /*! \fn
              * A function in order to update the gaps in amino acid chains of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param gaps The list of missing residues
              */
            void UpdateGapsInAminoAcidChains(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files, PdbPreprocessorMissingResidueVector gaps);
            /*! \fn
              * A function in order to update the gaps in amino acid chains of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param gaps The list of missing residues
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void UpdateGapsInAminoAcidChainsWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files,
                                                                    PdbPreprocessorMissingResidueVector gaps, int model_number = 1);
            /*! \fn
              * A function in order to access to the library residue by name from multiple library files
              * @param residue_name The name of a residue
              * @param lib_files The list of paths to the library files
              * @return library_residue
              */
/** @}*/
            /** \addtogroup Molecular_Data_Structure
               * @{
               */
            LibraryFileSpace::LibraryFileResidue* GetLibraryResidueByNameFromMultipleLibraryFiles(std::string residue_name, std::vector<std::string> lib_files);
            /*! \fn
              * A function in order to extract the alternate residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @return bool value
              */
/** @}*/
            /** \addtogroup Manipulators
               * @{
               */
            bool ExtractAlternateResidue(std::string pdb_file_path);
            /*! \fn
              * A function in order to extract the alternate residues of a pdb file
              * @param pdb_file The object of a pdb file
              * @return bool value
              */
            bool ExtractAlternateResidue(PdbFileSpace::PdbFile* pdb_file);
            /*! \fn
              * A function in order to remove unselected alternate residues of a pdb file
              * @param pdb_file The object of a pdb file
              * @param alternate_residue_map
              */
            void RemoveUnselectedAlternateResidues(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorAlternateResidueMap alternate_residue_map);
            /*! \fn
              * A function in order to remove unselected alternate residues of a pdb file
              * @param pdb_file The object of a pdb file
              * @param alternate_residue_map
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void RemoveUnselectedAlternateResiduesWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorAlternateResidueMap alternate_residue_map/*, int model_number = 1*/);
            /*! \fn
              * A function to do all preprocessing sequentially
              * @param pdb_file Pdb file object that has to be preprocessed
              * @param lib_files_path Paths of library files as database in order for preprocessing of the given pdb file
              * @param prep_files_path Paths of prep files as database in order for preprocessing of the given pdb file
              */
            void Preprocess(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files_path, std::vector<std::string> glycam_lib_files_path, std::vector<std::string> other_lib_files_path, std::vector<std::string> prep_files_path);
            /*! \fn
              * A function to apply all the updates on a pdb file
              * @param pdb_file A pdb file object that has to modified to reflect the updates
              * @param lib_files_path Paths of library files as database in order for preprocessing of the given pdb file
              */
            void ApplyPreprocessing(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files_path,
                                    std::vector<std::string> glycam_lib_files_path, std::vector<std::string> prep_files_path);
            /*! \fn
              * A function to apply all the updates on a pdb file
              * @param pdb_file A pdb file object that has to modified to reflect the updates
              * @param lib_files_path Paths of library files as database in order for preprocessing of the given pdb file
              */
            void ApplyPreprocessingWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files_path,
                                                           std::vector<std::string> glycam_lib_files_path, std::vector<std::string> prep_files_path, int model_number = 1);
            /*! \fn
              * A function to delete all to be deleted atoms and residues in a pdb file
              * @param pdb_file A pdb file object that has to modified to reflect the updates
              */
            void DeleteAllToBeDeletedEntities(PdbFileSpace::PdbFile* pdb_file);
            /*! \fn
              * A function to delete all to be deleted atoms and residues in a pdb file
              * @param pdb_file A pdb file object that has to modified to reflect the updates
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void DeleteAllToBeDeletedEntitiesWithTheGivenModelNumber(PdbFileSpace::PdbFile* pdb_file, int model_number = 1);
            /*! \fn
              * A function in order to extract the residues info of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @return bool value
              */
            bool ExtractResidueInfo(std::string pdb_file_path, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_filess, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to extract the residues info of a pdb file
              * @param pdb_file The object of a pdb file
              * @return bool value
              */
            bool ExtractResidueInfo(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function to calculate the overall charge of the model
              * @param pdb_file_path The path to the pdb file
              * @param lib_files Paths of library files as database in order for preprocessing of the given pdb file
              * @return model_charge Overal charge of the model
              */
            double CalculateModelCharge(std::string pdb_file_path, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function to calculate the overall charge of the model
              * @param pdb_file The object of a pdb file
              * @param lib_files Paths of library files as database in order for preprocessing of the given pdb file
              * @return model_charge Overal charge of the model
              */
            double CalculateModelCharge(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> amino_lib_files, std::vector<std::string> glycam_lib_files, std::vector<std::string> other_lib_files, std::vector<std::string> prep_files);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            PdbPreprocessorDisulfideBondVector disulfide_bonds_;                    /*!< List of detected disulfide bonds in a pdb file >*/
            PdbPreprocessorChainTerminationVector chain_terminations_;              /*!< List of chains detected in a pdb file >*/
            PdbPreprocessorHistidineMappingVector histidine_mappings_;              /*!< List of histidine residues detected in a pdb file >*/
            PdbPreprocessorMissingResidueVector missing_residues_;                  /*!< List of gaps (missing residues) detected in a pdb file >*/
            PdbPreprocessorUnrecognizedResidueVector unrecognized_residues_;        /*!< List of unrecognized residues detected in a pdb file >*/
            PdbPreprocessorRecognizedResidueVector recognized_residues_;            /*!< List of recognized residues detected in a pdb file >*/
            PdbPreprocessorUnrecognizedHeavyAtomVector unrecognized_heavy_atoms_;   /*!< List of unrecognized heavy atoms detected in a pdb file >*/
            PdbPreprocessorReplacedHydrogenVector replaced_hydrogens_;              /*!< List of removed/replaced hydrogen atoms detected in a pdb file >*/
            PdbPreprocessorAlternateResidueMap alternate_residue_map_;              /*!< Map of alternate residues detected in a pdb file >*/
            PdbPreprocessorToBeDeletedAtomVector to_be_deleted_atoms_;              /*!< List of atoms to be deleted in a pdb file >*/
            PdbPreprocessorToBeDeletedResidueVector to_be_deleted_residues_;        /*!< List of residues to be deleted in a pdb file >*/
            PdbPreprocessorResidueInfoMap residue_info_map_;                        /*!< Map of residues in a pdb file >*/

    };
}

#endif // PDBPREPROCESSOR_HPP
