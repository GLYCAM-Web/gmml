#ifndef PDBPREPROCESSOR_HPP
#define PDBPREPROCESSOR_HPP

#include <string>
#include <iostream>
#include <vector>

#include "../../../includes/FileSet/PdbFileSpace/pdbresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbfile.hpp"

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorDisulfideBond;
    class PdbPreprocessorChainTermination;
    class PdbPreprocessorHistidineMapping;
    class PdbPreprocessorMissingResidue;
    class PdbPreprocessorUnrecognizedResidue;
    class PdbPreprocessorUnrecognizedHeavyAtom;
    class PdbPreprocessorReplacedHydrogen;
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
            typedef std::vector<PdbPreprocessorUnrecognizedHeavyAtom*> PdbPreprocessorUnrecognizedHeavyAtomVector;
            typedef std::vector<PdbPreprocessorReplacedHydrogen*> PdbPreprocessorReplacedHydrogenVector;
            typedef std::vector<PdbFileSpace::PdbResidue*> PdbResidueVector;

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
              * An accessor function in order to access to the unrecognized heavy atoms
              * @return unrecognized_heavy_atoms_ attribute of the current object of this class
              */
            PdbPreprocessorUnrecognizedHeavyAtomVector GetUnrecognizedHeavyAtoms();
            /*! \fn
              * An accessor function in order to access to the replaced hydrogens
              * @return replaced_hydrogens_ attribute of the current object of this class
              */
            PdbPreprocessorReplacedHydrogenVector GetReplacedHydrogens();

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
            std::vector<std::string> GetUnrecognizedResidueNames(std::vector<std::string> pdb_residue_names, std::vector<std::string> dataset_residue_names);
            std::vector<std::string> GetRecognizedResidueNames(std::vector<std::string> pdb_residue_names, std::vector<std::string> dataset_residue_names);
            PdbResidueVector GetUnrecognizedResidues(PdbResidueVector pdb_residues, std::vector<std::string> unrecognized_residue_names);
            PdbResidueVector GetRecognizedResidues(PdbResidueVector pdb_residues, std::vector<std::string> recognized_residue_names);
            void ExtractUnrecognizedResidues(std::string pdb_file_path, std::vector<std::string> lib_files, std::vector<std::string> prep_files);
            std::vector<std::string> GetAllResidueNamesFromMultipleLibFiles(std::vector<std::string> lib_files);
            std::vector<std::string> GetAllResidueNamesFromMultiplePrepFiles(std::vector<std::string> prep_files);
            std::vector<std::string> GetAllResidueNamesFromDatasetFiles(std::vector<std::string> lib_files, std::vector<std::string> prep_files);
            PdbResidueVector GetAllCYSResidues(PdbResidueVector pdb_residues);
            double GetDistanceofCYS(PdbFileSpace::PdbResidue* first_residue, PdbFileSpace::PdbResidue* second_residue, PdbFileSpace::PdbFile* pdb_file, PdbFileSpace::PdbFile::PdbResidueAtomsMap residue_atom_map);
            void ExtractCYSResidues(std::string pdb_file_path);
            PdbResidueVector GetAllHISResidues(PdbResidueVector pdb_residues);
            void ExtractHISResidues(std::string pdb_file_path);

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
            PdbPreprocessorUnrecognizedHeavyAtomVector unrecognized_heavy_atoms_;
            PdbPreprocessorReplacedHydrogenVector replaced_hydrogens_;

    };
}

#endif // PDBPREPROCESSOR_HPP
