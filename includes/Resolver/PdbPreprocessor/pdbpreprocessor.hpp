#ifndef PDBPREPROCESSOR_HPP
#define PDBPREPROCESSOR_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../../../includes/FileSet/PdbFileSpace/pdbresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"

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
    class PdbPreprocessor
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<PdbPreprocessorDisulfideBond*> PdbPreprocessorDisulfideBondVector;
            typedef std::vector<PdbPreprocessorChainTermination*> PdbPreprocessorChainTerminationVector;
            typedef std::vector<PdbPreprocessorHistidineMapping*> PdbPreprocessorHistidineMappingVector;
            typedef std::vector<PdbPreprocessorMissingResidue*> PdbPreprocessorMissingResidueVector;
            typedef std::vector<PdbPreprocessorUnrecognizedResidue*> PdbPreprocessorUnrecognizedResidueVector;
            typedef std::vector<PdbPreprocessorUnrecognizedResidue*> PdbPreprocessorRecognizedResidueVector;
            typedef std::vector<PdbPreprocessorUnrecognizedHeavyAtom*> PdbPreprocessorUnrecognizedHeavyAtomVector;
            typedef std::vector<PdbPreprocessorReplacedHydrogen*> PdbPreprocessorReplacedHydrogenVector;
            typedef std::map<char, std::vector<int> > PdbPreprocessorChainIdSequenceNumbersMap;
            typedef std::map<char, std::vector<char> > PdbPreprocessorChainIdInsertionCodeMap;
            typedef std::map<std::string, PdbPreprocessorAlternateResidue*> PdbPreprocessorAlternateResidueMap;


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

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function in order to access to the unrecognized residue names of pdb
              * @param pdb_residue_names The list of residue names in the current pdb file
              * @param data_set_residue_names The list of residue names from library and prep files
              * @return unrecognized_residue_names
              */
            std::vector<std::string> GetUnrecognizedResidueNames(std::vector<std::string> pdb_residue_names, std::vector<std::string> dataset_residue_names);
            /*! \fn
              * A function in order to access to the recognized residue names of pdb
              * @param pdb_residue_names The list of residue names in the current pdb file
              * @param data_set_residue_names The list of residue names from library and prep files
              * @return recognized_residue_names
              */
            std::vector<std::string> GetRecognizedResidueNames(std::vector<std::string> pdb_residue_names, std::vector<std::string> dataset_residue_names);
            /*! \fn
              * A function in order to access to the unrecognized residues of pdb
              * @param pdb_residues The list of residues in the current pdb file
              * @param unrecognized_residue_names The list of unrecognized residue names
              * @return unrecognized_residues
              */
            PdbFileSpace::PdbFile::PdbResidueVector GetUnrecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, std::vector<std::string> unrecognized_residue_names);
            /*! \fn
              * A function in order to access to the recognized residues of pdb
              * @param pdb_residues The list of residues in the current pdb file
              * @param recognized_residue_names The list of recognized residue names
              * @return recognized_residues
              */
            PdbFileSpace::PdbFile::PdbResidueVector GetRecognizedResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues, std::vector<std::string> recognized_residue_names);

            /*! \fn
              * A function in order to access to the list of all residue names from lib files
              * @param lib_files The list of paths to library files
              * @return all_residue_names
              */
            std::vector<std::string> GetAllResidueNamesFromMultipleLibFiles(std::vector<std::string> lib_files);
            /*! \fn
              * A function in order to access to the list of all residue names from prep files
              * @param prep_files The list of paths to prep files
              * @return all_residue_names
              */
            std::vector<std::string> GetAllResidueNamesFromMultiplePrepFiles(std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to access to the list of all residue names from library and prep files
              * @param lib_files The list of paths to library files
              * @param prep_files The list of paths to prep files
              * @return all_residue_names
              */
            std::vector<std::string> GetAllResidueNamesFromDatasetFiles(std::vector<std::string> lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to extract the unrecognized residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              */
            void ExtractUnrecognizedResidues(std::string pdb_file_path, std::vector<std::string> lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to remove the unrecognized residues of a pdb file
              * @param pdb_file The object of a pdb file
              * @param unknown_heavy_atoms The list of unrecognized residues
              */
            void RemoveUnrecognizedResidues(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorUnrecognizedResidueVector unrecognized_residues);

            /*! \fn
              * A function in order to extract the recognized residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              */
            void ExtractRecognizedResidues(std::string pdb_file_path, std::vector<std::string> lib_files, std::vector<std::string> prep_files);

            /*! \fn
              * A function in order to access to the list of CYS residues
              * @param pdb_residues The list of pdb residues
              * @return all_cys_residues
              */
            PdbFileSpace::PdbFile::PdbResidueVector GetAllCYSResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues);
            /*! \fn
              * A function in order to access to the distance of a pair of CYS residues
              * @param first_residue The first residue of CYS pair
              * @param second_residue The second residue of CYS pair
              * @param pdb_file Pdb file object
              * @param residue_atom_map A map between a residue and its atoms
              * @return distance
              */
            double GetDistanceofCYS(PdbFileSpace::PdbResidue* first_residue, PdbFileSpace::PdbResidue* second_residue, PdbFileSpace::PdbFile* pdb_file, PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atom_map);
            /*! \fn
              * A function in order to extract the CYS residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              */
            void ExtractCYSResidues(std::string pdb_file_path);
            /*! \fn
              * A function in order to remove the unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param disulfide_bonds The list of disulfide bonds
              */
            void UpdateCYSResidues(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorDisulfideBondVector disulfide_bonds);

            /*! \fn
              * A function in order to access to the list of HIS residues
              * @param pdb_residues The list of pdb residues
              * @return all_his_residues
              */
            PdbFileSpace::PdbFile::PdbResidueVector GetAllHISResidues(PdbFileSpace::PdbFile::PdbResidueVector pdb_residues);
            /*! \fn
              * A function in order to extract the HIS residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              */
            void ExtractHISResidues(std::string pdb_file_path);
            /*! \fn
              * A function in order to update histidine mapping of a pdb file
              * @param pdb_file The object of a pdb file
              * @param histidine_mappings The list of histidine mappings
              */
            void UpdateHISMapping(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorHistidineMappingVector histidine_mappings);

            /*! \fn
              * A function in order to access to the list of unknown heavy atoms of a residue
              * @param pdb_atom_names_of_residue The list of atom names of a pdb residue
              * @param dataset_atom_names_of_residue The list of atom names of a residue from dataset
              * @return unknown_heavy_atom_names_of_residue
              */
            std::vector<std::string> GetUnknownHeavyAtomNamesOfResidue(std::vector<std::string> pdb_atom_names_of_residue, std::vector<std::string> dataset_atom_names_of_residue);
            /*! \fn
              * A function in order to access to the list of atom names of a residue from multiple library files
              * @param residue_name The name of the residue
              * @param lib_files The list of paths to the library files
              * @return all_atom_names_of_residue
              */
            std::vector<std::string> GetAllAtomNamesOfResidueFromMultipleLibFiles(std::string residue_name, std::vector<std::string> lib_files);
            /*! \fn
              * A function in order to access to the list of atom names of a residue from multiple prep files
              * @param residue_name The name of the residue
              * @param prep_files The list of paths to the prep files
              * @return all_atom_names_of_residue
              */
            std::vector<std::string> GetAllAtomNamesOfResidueFromMultiplePrepFiles(std::string residue_name, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to access to the list of all atom names of a residue from library and prep files
              * @param lib_files The list of paths to library files
              * @param prep_files The list of paths to prep files
              * @return all_atom_names
              */
            std::vector<std::string> GetAllAtomNamesOfResidueFromDatasetFiles(std::string residue_name, std::vector<std::string> lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to access to the unknown heavy atoms of a residue
              * @param pdb_atoms The list of pbd atoms
              * @param dataset_atom_names_of_residue The list of atom names of a residue from dataset
              * @return unknown_heavy_atoms_of_residue
              */
            PdbFileSpace::PdbFile::PdbAtomVector GetUnknownHeavyAtomsOfResidue(PdbFileSpace::PdbFile::PdbAtomVector pdb_atoms, std::vector<std::string> dataset_atom_names_of_residue);
            /*! \fn
              * A function in order to extract the unknown heavy atoms of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              */
            void ExtractUnknownHeavyAtoms(std::string pdb_file_path, std::vector<std::string> lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to remove the unknown heavy atoms of a pdb file
              * @param pdb_file The object of a pdb file
              * @param unknown_heavy_atoms The list of unknown heavy atoms
              */
            void RemoveUnknownHeavyAtoms(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms);

            /*! \fn
              * A function in order to access to the removed hydrogen names of a residue
              * @param pdb_atom_names_of_residue The list of atom names of a residue from pdb file
              * @param dataset_atom_names_of_residue The list of atom names of a residue from dataset
              * @return removed_hydrogen_names_of_residue
              */
            std::vector<std::string> GetRemovedHydrogenNamesOfResidue(std::vector<std::string> pdb_atom_names_of_residue, std::vector<std::string> dataset_atom_names_of_residue);
            /*! \fn
              * A function in order to access to the removed hydrogens of a residue
              * @param pdb_atoms The list of pdb atoms
              * @param dataset_atom_names_of_residue The list of atom names of a residue from dataset
              * @return removed_hydrogens_of_residue
              */
            PdbFileSpace::PdbFile::PdbAtomVector GetRemovedHydrogensOfResidue(PdbFileSpace::PdbFile::PdbAtomVector pdb_atoms, std::vector<std::string> dataset_atom_names_of_residue);
            /*! \fn
              * A function in order to extract the removed hydrogens of a pdb file
              * @param pdb_file_path The path to the pdb file
              * @param lib_files The list of paths to the library files
              * @param prep_files The list of paths to the prep files
              */
            void ExtractRemovedHydrogens(std::string pdb_file_path, std::vector<std::string> lib_files, std::vector<std::string> prep_files);
            /*! \fn
              * A function in order to remove the removed hydrogens of a pdb file
              * @param pdb_file The object of a pdb file
              * @param replaced_hydrogens The list of replaced hydrogens
              */
            void RemoveRemovedHydrogens(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorReplacedHydrogenVector replaced_hydrogens);

            /*! \fn
              * A function in order to extract the amino acid chains of a pdb file
              * @param pdb_file_path The path to the pdb file
              */
            void ExtractAminoAcidChains(std::string pdb_file_path);
            /*! \fn
              * A function in order to update the amino acid chains of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param gaps The list of chain terminations
              */
            void UpdateAminoAcidChains(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> lib_files, PdbPreprocessorChainTerminationVector chain_termination);

            /*! \fn
              * A function in order to extract the gaps in amino acid chains of a pdb file
              * @param pdb_file_path The path to the pdb file
              */
            void ExtractGapsInAminoAcidChains(std::string pdb_file_path);
            /*! \fn
              * A function in order to update the gaps in amino acid chains of a pdb file
              * @param pdb_file The object of a pdb file
              * @param lib_files The list of paths to the library files
              * @param gaps The list of missing residues
              */
            void UpdateGapsInAminoAcidChains(PdbFileSpace::PdbFile* pdb_file, std::vector<std::string> lib_files, PdbPreprocessorMissingResidueVector gaps);

            /*! \fn
              * A function in order to access to the library residue by name from multiple library files
              * @param residue_name The name of a residue
              * @param lib_files The list of paths to the library files
              * @return library_residue
              */
            LibraryFileSpace::LibraryFileResidue* GetLibraryResidueByNameFromMultipleLibraryFiles(std::string residue_name, std::vector<std::string> lib_files);

            /*! \fn
              * A function in order to extract the alternate residues of a pdb file
              * @param pdb_file_path The path to the pdb file
              */
            void ExtractAlternateResidue(std::string pdb_file_path);
            /*! \fn
              * A function in order to remove unselected alternate residues of a pdb file
              * @param pdb_file The object of a pdb file
              * @param alternate_residue_map
              */
            void RemoveUnselectedAlternateResidues(PdbFileSpace::PdbFile* pdb_file, PdbPreprocessorAlternateResidueMap alternate_residue_map);

            void RetrievePreprocessingInformation(std::string pdb_file_path, std::vector<std::string> lib_files_path, std::vector<std::string> prep_files_path);



            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            PdbPreprocessorDisulfideBondVector disulfide_bonds_;
            PdbPreprocessorChainTerminationVector chain_terminations_;
            PdbPreprocessorHistidineMappingVector histidine_mappings_;
            PdbPreprocessorMissingResidueVector missing_residues_;
            PdbPreprocessorUnrecognizedResidueVector unrecognized_residues_;
            PdbPreprocessorRecognizedResidueVector recognized_residues_;
            PdbPreprocessorUnrecognizedHeavyAtomVector unrecognized_heavy_atoms_;
            PdbPreprocessorReplacedHydrogenVector replaced_hydrogens_;
            PdbPreprocessorAlternateResidueMap alternate_residue_map_;

    };
}

#endif // PDBPREPROCESSOR_HPP
