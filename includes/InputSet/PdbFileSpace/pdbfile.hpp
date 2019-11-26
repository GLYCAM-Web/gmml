//Author: Alireza Khatamian
//Modified by: Dave Montgomery

#ifndef PDBFILE_HPP
#define PDBFILE_HPP

#include "inputfile.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace PdbFileSpace
{
    //Title Section
    class PdbHeaderCard;
    class PdbObsoleteSection;
    class PdbTitleSection;
    class PdbSplitSection;
    class PdbCaveatSection;
    class PdbCompoundSection;
    class PdbSourceSection;
    class PdbSourceCard;
    class PdbKeywordsSection;
    class PdbExperimentalDataSection;
    class PdbNumModelCard;
    class PdbModelTypeSection;
    class PdbAuthorSection;
    class PdbRevisionDataSection;
    class PdbRevisionDataCard;
    class PdbSupersededEntriesSection;
    class PdbJournalSection;
    class PdbRemarkSection;

    //Primary Structure Section
    class PdbDatabaseReferenceSection;
    class PdbDatabaseReference;
    class PdbSequenceAdvancedSection;
    class PdbSequenceAdvancedCard;
    class PdbResidueSequenceSection;
    class PdbResidueModificationSection;
    class PdbResidue;

    //Heterogen Section
    class PdbHeterogenSection;
    class PdbHeterogenNameSection;
    class PdbHeterogenSynonymSection;
    class PdbFormulaSection;

    //Secondary Structure Section
    class PdbHelixSection;
    class PdbSheetSection;

    //Connectivity Annotation Section
    class PdbDisulfideBondSection;
    class PdbLinkCard;
    class PdbLinkSection;
    class PdbCISPeptideSection;
    
    class PdbCISPeptideCard;

    //Miscellaneous Features Section
    class PdbSiteSection;

    //Crystallographic and Coordinate Transformation Section
    class PdbCrystallographicCard;
    class PdbOriginXnSection;
    class PdbScaleNSection;
    class PdbMatrixNSection;

    //Coordinate Section
    class PdbModelSection;
    class PdbModelCard;
    class PdbAtomSection;
    class PdbAtomCard;
    //class PdbAnisotropicSection;
    //class PdbTerminalSection;

    //Connectivity Section
    class PdbConnectSection;

    //Bookkeeping Section
    class PdbMasterCard;
    // class PdbEndCard;
    


    class PdbFile: public InputFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of sources
              */
            typedef std::vector<PdbSourceCard*> SourceCardVector;
            /*! \typedef
              * List of sequence advanced cards
              */
            typedef std::vector<PdbSequenceAdvancedCard*> SequenceAdvancedCardVector;
            /*! \typedef
              * List of revision_datas
              */
            typedef std::vector<PdbRevisionDataCard*> RevisionDataCardVector;
            /*! \typedef
              * List of CIS peptide cards
              */
            typedef std::vector<PdbCISPeptideCard*> CISPeptideCardVector;
            /*! \typedef
              * List of database references
              */
            typedef std::vector<PdbDatabaseReference*> DatabaseReferenceVector;
            /*! \typedef
              * List of residues
              */
            typedef std::vector<PdbResidue*> PdbResidueVector;
            /*! \typedef
              * List of pdb atom
              */
            typedef std::vector<PdbAtomCard*> PdbAtomCardVector;
            /*! \typedef
              * A mapping between a
              */
            typedef std::map<std::string, PdbAtomCardVector* > PdbResidueAtomsMap;
            /*! \typedef
              * Mapping between old serial number and new one that has been changed during a process
              */
            typedef std::map<int, int> PdbSerialNumberMapping;
            /*! \typedef
              * Mapping between old sequence number and new one that has been changed during a process
              */
            typedef std::map<int, int> PdbSequenceNumberMapping;

            typedef std::vector<std::pair<char, int> > PdbPairVectorTerCardPositions;
            typedef std::vector<std::pair<std::string, std::string> > PdbPairVectorAtomNamePositionFlag;
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbFile();
            /*! \fn
              * Constructor
              * @param pdb_file An existing pdb file path to be read
              */
            PdbFile(const std::string& pdb_file);
            /*! \fn
              * Constructor
              * @param atomStream A stringstream of atom cards
              */
            PdbFile(std::stringstream& atomStream);
            
            /*! \fn
              * Load PDB file
              */
            PdbFile* LoadPdbFile();
            /*! \fn
              * @param pdb_file An existing pdb file path to be read
              */
            PdbFile* LoadPdbFile(const std::string& pdb_file);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to pdb file path of the current object
              * @return path_ attribute of the current object of this class
              */
            std::string GetPath();
            /*! \fn
              * An accessor function in order to access to the header of the current object
              * @return header_ attribute of the current object of this class
              */
            PdbHeaderCard* GetHeader();
            /*! \fn
              * An accessor function in order to access to the obsolete of the current object
              * @return obsolete_ attribute of the current object of this class
              */
            PdbObsoleteSection* GetObsoleteCards();
            /*! \fn
              * An accessor function in order to access to the title of the current object
              * @return title_ attribute of the current object of this class
              */
            PdbTitleSection* GetTitle();
            /*! \fn
              * An accessor function in order to access to the split of the current object
              * @return split_ attribute of the current object of this class
              */
            PdbSplitSection* GetSplit();
            /*! \fn
              * An accessor function in order to access to the caveat of the current object
              * @return caveat_ attribute of the current object of this class
              */
            PdbCaveatSection* GetCaveat();
            /*! \fn
              * An accessor function in order to access to the compound attribute of the current object
              * @return compound_ attribute of the current object of this class
              */
            PdbCompoundSection* GetCompound();
            /*! \fn
              * An accessor function in order to access to the source of the current object
              * @return source_ attribute of the current object of this class
              */
            PdbSourceSection* GetSourceCards();
            /*! \fn
              * An accessor function in order to access to the keywords of the current object
              * @return keywords_ attribute of the current object of this class
              */
            PdbKeywordsSection* GetKeywords();
            /*! \fn
              * An accessor function in order to access to the experimental_data of the current object
              * @return experimental_data_ attribute of the current object of this class
              */
            PdbExperimentalDataSection* GetExperimentalData();
            /*! \fn
              * An accessor function in order to access to the number of models attribute of the current object
              * @return number_of_models_ attribute of the current object of this class
              */
            PdbNumModelCard* GetNumberOfModels();
            /*! \fn
              * An accessor function in order to access to the model type attribute of the current object
              * @return model_type_ attribute of the current object of this class
              */
            PdbModelTypeSection* GetModelType();
            /*! \fn
              * An accessor function in order to access to the author of the current object
              * @return author_ attribute of the current object of this class
              */
            PdbAuthorSection* GetAuthor();
            /*! \fn
              * An accessor function in order to access to the revision_data of the current object
              * @return revision_data_ attribute of the current object of this class
              */
            PdbRevisionDataSection* GetRevisionDataCards();
            /*! \fn
              * An accessor function in order to access to the superseded_entries of the current object
              * @return superseded_entries_ attribute of the current object of this class
              */
            PdbSupersededEntriesSection* GetSupersededEntriesCards();
            /*! \fn
              * An accessor function in order to access to the journal of the current object
              * @return journal_ attribute of the current object of this class
              */
            PdbJournalSection* GetJournal();
            /*! \fn
              * An accessor function in order to access to the remark_cards type attribute of the current object
              * @return remark_cards_ attribute of the current object of this class
              */
            PdbRemarkSection* GetRemarks();
            /*! \fn
              * An accessor function in order to access to the database_reference of the current object
              * @return database_reference_ attribute of the current object of this class
              */
            PdbDatabaseReferenceSection* GetDatabaseReferences();
            /*! \fn
              * An accessor function in order to access to the sequence_advanced of the current object
              * @return sequence_advanced_ attribute of the current object of this class
              */
            PdbSequenceAdvancedSection* GetSequenceAdvanced();
            /*! \fn
              * An accessor function in order to access to the residue sequence attribute of the current object
              * @return residue_sequence_ attribute of the current object of this class
              */
            PdbResidueSequenceSection* GetResiduesSequence();
            /*! \fn
              * An accessor function in order to access to the residue modification attribute of the current object
              * @return residue_modification_cards_ attribute of the current object of this class
              */
            PdbResidueModificationSection* GetResidueModification();
            /*! \fn
              * An accessor function in order to access to the heterogens attribute of the current object
              * @return heterogen_cards_ attribute of the current object of this class
              */
            PdbHeterogenSection* GetHeterogenCards();
            /*! \fn
              * An accessor function in order to access to the heterogens name attribute of the current object
              * @return heterogens_name_ attribute of the current object of this class
              */
            PdbHeterogenNameSection* GetHeterogenNameCards();
            /*! \fn
              * An accessor function in order to access to the heterogen synonyms attribute of the current object
              * @return heterogen_synonyms_ attribute of the current object of this class
              */
            PdbHeterogenSynonymSection* GetHeterogenSynonymCards();
            /*! \fn
              * An accessor function in order to access to the formulas attribute of the current object
              * @return formulas_ attribute of the current object of this class
              */
            PdbFormulaSection* GetFormulaCards();
            /*! \fn
              * An accessor function in order to access to the helixes attribute of the current object
              * @return helix_cards_ attribute of the current object of this class
              */
            PdbHelixSection* GetHelixCards();
            /*! \fn
              * An accessor function in order to access to the sheet_cards attribute of the current object
              * @return sheet_cards_ attribute of the current object of this class
              */
            PdbSheetSection* GetSheets();
            /*! \fn
              * An accessor function in order to access to the disulfide bonds attribute of the current object
              * @return disulfide_bonds_ attribute of the current object of this class
              */
            PdbDisulfideBondSection* GetDisulfideBonds();
            /*! \fn
              * An accessor function in order to access to the links attribute of the current object
              * @return link_cards_ attribute of the current object of this class
              */
            PdbLinkSection* GetResidueLinkCards();
            /*! \fn
              * An accessor function in order to access to the cis_peptide of the current object
              * @return cis_peptide_ attribute of the current object of this class
              */
            PdbCISPeptideSection* GetCISPeptide();
            /*! \fn
              * An accessor function in order to access to the site_cards attribute of the current object
              * @return site_cards_ attribute of the current object of this class
              */
            PdbSiteSection* GetSites();
            /*! \fn
              * An accessor function in order to access to the crystallography attribute of the current object
              * @return crystallography_ attribute of the current object of this class
              */
            PdbCrystallographicCard* GetCrystallography();
            /*! \fn
              * An accessor function in order to access to the origins attribute of the current object
              * @return origins_ attribute of the current object of this class
              */
            PdbOriginXnSection* GetOrigins();
            /*! \fn
              * An accessor function in order to access to the scales attribute of the current object
              * @return scales_ attribute of the current object of this class
              */
            PdbScaleNSection* GetScales();
            /*! \fn
              * An accessor function in order to access to the matrices attribute of the current object
              * @return matrices_ attribute of the current object of this class
              */
            PdbMatrixNSection* GetMatrices();
            /*! \fn
              * An accessor function in order to access to the models attribute of the current object
              * @return models_ attribute of the current object of this class
              */
            PdbModelSection* GetModels();
            /*! \fn
              * An accessor function in order to access to the connectivities attribute of the current object
              * @return connectivities_ attribute of the current object of this class
              */
            PdbConnectSection* GetConnectivities();
            /*! \fn
              * An accessor function in order to access to the serial number mapping attribute of the current object
              * @return serial_number_mapping attribute of the current object of this class
              */
            PdbSerialNumberMapping GetSerialNumberMapping();
            /*! \fn
              * An accessor function in order to access to the sequence number mapping attribute of the current object
              * @return sequence_number_mapping attribute of the current object of this class
              */
            PdbSequenceNumberMapping GetSequenceNumberMapping();
            /*! \fn
              * An accessor function in order to access to all residue names of the current object
              * @return residue_names All residue names of the current object of this class
              */
            PdbPairVectorAtomNamePositionFlag GetAllResidueNames();
            /*! \fn
              * An accessor function in order to access to all residue names from atom card of the current object
              * @return residue_names All residue names from atom card of the current object of this class
              */
            PdbPairVectorAtomNamePositionFlag GetAllResidueNamesFromAtomSection();
            /*! \fn
              * An accessor function in order to access to all residues of the current object
              * @return all_residues_ All residues of the current object of this class
              */
            PdbResidueVector GetAllResidues();
            /*! \fn
              * An accessor function in order to access to all residues from atom card of the current object
              * @return residues All resdidues from atom card of the current object of this class
              */
            PdbResidueVector GetAllResiduesFromAtomSection();
            /*! \fn
              * An accessor function in order to access to all atoms of a residue of the current object
              * @param residue The given residue to return all the atoms of it
              * @return atoms_of_residue All atoms of a resdidue of the current object of this class
              */
            PdbAtomCardVector GetAllAtomsOfResidue(PdbResidue* residue);
            /*! \fn
              * An accessor function in order to access to all atoms of all residues of the current object
              * @return residue_atom_map The map between all residues and their atoms of the current object of this class
              */
            PdbResidueAtomsMap GetAllAtomsOfResidues();
            /*! \fn
              * An accessor function in order to access to all atoms of all residues of the current object
              * @return atoms The vector of all atoms in a pdb file in the order that they appear in the file
              */
            PdbResidueAtomsMap GetAllAtomsInOrder(std::vector<std::string>& key_order);
            /*! \fn
              * An accessor function in order to access to atom of a residue of the current object using the atom name
              * @param residue The given residue to return one of its atoms
              * @param atom_name The atom name of the desired atom object
              * @param residue_atom_map The map between residues and their atoms
              * @return atom The atom object of the given residue of the current object of this class
              */
            PdbAtomCard* GetAtomOfResidueByName(PdbResidue* residue, std::string atom_name, PdbResidueAtomsMap residue_atom_map);
            /*! \fn
              * An accessor function in order to access to atom of a residue of the current object using the atom name
              * @param residue The given residue to return one of its atoms
              * @param atom_name The atom name of the desired atom object
              * @return atom The atom object of the given residue of the current object of this class
              */
            PdbAtomCard* GetAtomOfResidueByName(PdbResidue* residue, std::string atom_name);
            /*! \fn
              * An accessor function in order to access to atom of a residue of the current object using the atom key
              * @param atom_key The atom key of the desired atom object
              * @return heterogen_atom_map The map between heterogen atoms of the current object of this class and their keys
              */
            PdbAtomCard* GetAtomOfResidueByAtomKey(std::string atom_key);
            /*! \fn
              * An accessor function in order to access to atom of a residue of the current object using the serial number
              * @param serial_number The serial number attribute of the desired atom object
              * @return heterogen_atom_map The map between heterogen atoms of the current object of this class and their keys
              */
            PdbAtomCard* GetAtomBySerialNumber(int serial_number);
            /*! \fn
              * An accessor function in order to access to all atom names of a residue of the current object
              * @param residue The given residue to return all of its atom names
              * @param residue_atom_map The map between residues and their atoms
              * @return atom_names The atom names of the given residue of the current object of this class
              */
            std::vector<std::string> GetAllAtomNamesOfResidue(PdbResidue* residue, PdbResidueAtomsMap residue_atom_map);
            /*! \fn
              * An accessor function in order to access to all atom names of a residue of the current object
              * @param residue The given residue to return all of its atom names
              * @return atom_names The atom names of the given resdidue of the current object of this class
              */
            std::vector<std::string> GetAllAtomNamesOfResidue(PdbResidue* residue);
            /*! \fn
              * A finction that returns chain id and sequence number of a residue that has been placed after a residue that has no tail or has at least 2 tails
              * @return A 2-D vector of chain id and sequence number of residues that have been placed after residues that have no tail or have at least 2 tails
              */
            PdbPairVectorTerCardPositions GetAllTerCardPositions(std::vector<std::string> glycam_residue_names);
            /*! \fn
              * An accessor function in order to access tthe master card of the current object
              * @return master_ The atom names of the given resdidue of the current object of this class
              */
            PdbMasterCard* GetMasterCard();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A function in order to set the path of the pdb file
              * @param pdb_path Pdb file path
              */
            void SetPath(std::string pdb_path);
            /*! \fn
              * A mutator function in order to set the header card of the current object
              * Set the header_ attribute of the current pdb file
              * @param header The header attribute of the current object
              */
            void SetHeader(PdbHeaderCard* header);
            /*! \fn
              * A mutator function in order to set the obsolete section of the current object
              * Set the obsolete_ attribute of the current pdb file
              * @param obsolete The obsolete attribute of the current object
              */
            void SetObsolete(PdbObsoleteSection* obsolete);
            /*! \fn
              * A mutator function in order to set the title card of the current object
              * Set the title_ attribute of the current pdb file
              * @param title The title attribute of the current object
              */
            void SetTitle(PdbTitleSection* title);
            /*! \fn
              * A mutator function in order to set the split card of the current object
              * Set the split_ attribute of the current pdb file
              * @param split The split attribute of the current object
              */
            void SetSplit(PdbSplitSection* split);
            /*! \fn
              * A mutator function in order to set the caveat card of the current object
              * Set the caveat_ attribute of the current pdb file
              * @param caveat The caveat attribute of the current object
              */
            void SetCaveat(PdbCaveatSection* caveat);
            /*! \fn
              * A mutator function in order to set the compound card of the current object
              * Set the compound_ attriSourcebute of the current pdb file
              * @param compound The compound attribute of the current object
              */
            void SetCompound(PdbCompoundSection* compound);
            /*! \fn
              * A mutator function in order to set the source card of the current object
              * Set the source_ attribute of the current pdb file
              * @param source The source attribute of the current object
              */
            void SetSourceCards(PdbSourceSection* source);
            /*! \fn
              * A mutator function in order to set the keywords card of the current object
              * Set the keywords_ attribute of the current pdb file
              * @param keywords The keywords attribute of the current object
              */
            void SetKeywords(PdbKeywordsSection* keywords);
            /*! \fn
              * A mutator function in order to set the experimental_data card of the current object
              * Set the experimental_data_ attribute of the current pdb file
              * @param experimental_data The experimental_data attribute of the current object
              */
            void SetExperimentalData(PdbExperimentalDataSection* experimental_data);
            /*! \fn
              * A mutator function in order to set the number of models card of the current object
              * Set the number_of_models_ attribute of the current pdb file
              * @param number_of_models The number of models attribute of the current object
              */
            void SetNumberOfModels(PdbNumModelCard* number_of_models);
            /*! \fn
              * A mutator function in order to set the model_type card of the current object
              * Set the model_type_ attribute of the current pdb file
              * @param model_type The model type attribute of the current object
              */
            void SetModelType(PdbModelTypeSection* model_type);
            /*! \fn
              * A mutator function in order to set the author card of the current object
              * Set the author_ attribute of the current pdb file
              * @param author The author attribute of the current object
              */
            void SetAuthor(PdbAuthorSection* author);
            /*! \fn
              * A mutator function in order to set the revision_data card of the current object
              * Set the revision_data_ attribute of the current pdb file
              * @param revision_data The revision_data attribute of the current object
              */
            void SetRevisionDataCards(PdbRevisionDataSection* revision_data);
            /*! \fn
              * A mutator function in order to set the superseded_entries card of the current object
              * Set the superseded_entries_ attribute of the current pdb file
              * @param superseded_entries The superseded_entries attribute of the current object
              */
            void SetSupersededEntriesCards(PdbSupersededEntriesSection* superseded_entries);
            /*! \fn
              * A mutator function in order to set the journal of the current object
              * Set the journal_ attribute of the current pdb file
              * @param journal The journal of the current object
              */
            void SetJournal(PdbJournalSection* journal);
            /*! \fn
              * A mutator function in order to set the remark card of the current object
              * Set the remark_cards_ attribute of the current pdb file
              * @param remark_cards The model type attribute of the current object
              */
            void SetRemarks(PdbRemarkSection* remark_cards);
            /*! \fn
              * A mutator function in order to set the database_reference card of the current object
              * Set the database_reference_ attribute of the current pdb file
              * @param database_reference The database_reference attribute of the current object
              */
            void SetDatabaseReferences(PdbDatabaseReferenceSection* database_reference);
            /*! \fn
              * A mutator function in order to set the sequence_advanced card of the current object
              * Set the sequence_advanced_ attribute of the current pdb file
              * @param sequence_advanced The sequence_advanced attribute of the current object
              */
            void SetSequenceAdvanced(PdbSequenceAdvancedSection* sequence_advanced);
            /*! \fn
              * A mutator function in order to set the residues sequence card of the current object
              * Set the residues_sequence_ attribute of the current pdb file
              * @param residues_sequence The residues sequence attribute of the current object
              */
            void SetResiduesSequence(PdbResidueSequenceSection* residues_sequence);
            /*! \fn
              * A mutator function in order to set the residue modification card of the current object
              * Set the residue_modification_cards_ attribute of the current pdb file
              * @param residue_modification_cards The residue modification attribute of the current object
              */
            void SetResidueModification(PdbResidueModificationSection* residue_modification_cards);
            /*! \fn
              * A mutator function in order to set the heterogens card of the current object
              * Set the heterogen_cards_ attribute of the current pdb file
              * @param heterogens The heterogens attribute of the current object
              */
            void SetHeterogens(PdbHeterogenSection* heterogen_cards);
            /*! \fn
              * A mutator function in order to set the heterogens name card of the current object
              * Set the heterogens_name_ attribute of the current pdb file
              * @param heterogens_name The heterogens name attribute of the current object
              */
            void SetHeterogensName(PdbHeterogenNameSection* heterogens_name);
            /*! \fn
              * A mutator function in order to set the heterogen synonyms card of the current object
              * Set the heterogen_synonyms_ attribute of the current pdb file
              * @param heterogen_synonym_cards The heterogen synonyms attribute of the current object
              */
            void SetHeterogenSynonyms(PdbHeterogenSynonymSection* heterogen_synonym_cards);
            /*! \fn
              * A mutator function in order to set the formulas card of the current object
              * Set the formulas_ attribute of the current pdb file
              * @param formulas The formulas attribute of the current object
              */
            void SetFormulas(PdbFormulaSection* formulas);
            /*! \fn
              * A mutator function in order to set the helixes card of the current object
              * Set the helix_cards_ attribute of the current pdb file
              * @param helixes The helixes attribute of the current object
              */
            void SetHelixes(PdbHelixSection* helixes);
            /*! \fn
              * A mutator function in order to set the sheet_cards card of the current object
              * Set the sheet_cards_ attribute of the current pdb file
              * @param sheet_cards The sheet_cards attribute of the current object
              */
            void SetSheets(PdbSheetSection* sheet_cards);
            /*! \fn
              * A mutator function in order to set the disulfide bonds card of the current object
              * Set the disulfide_bonds_ attribute of the current pdb file
              * @param disulfide_bonds The disulfide bonds attribute of the current object
              */
            void SetDisulfideBonds(PdbDisulfideBondSection* disulfide_bonds);
            /*! \fn
              * A mutator function in order to set the links card of the current object
              * Set the link_cards_ attribute of the current pdb file
              * @param links The links attribute of the current object
              */
            void SetLinks(PdbLinkSection* links);
            /*! \fn
              * A mutator function in order to set the cis_peptide card of the current object
              * Set the cis_peptide_ attribute of the current pdb file
              * @param cis_peptide The cis_peptide attribute of the current object
              */
            void SetCISPeptide(PdbCISPeptideSection* cis_peptide);
            /*! \fn
              * A mutator function in order to set the site_cards card of the current object
              * Set the site_cards_ attribute of the current pdb file
              * @param site_cards The compound attribute of the current object
              */
            void SetSites(PdbSiteSection* site_cards);
            /*! \fn
              * A mutator function in order to set the crystallography card of the current object
              * Set the crystallography_ attribute of the current pdb file
              * @param crystallography The crystallography attribute of the current object
              */
            void SetCrystallography(PdbCrystallographicCard* crystallography);
            /*! \fn
              * A mutator function in order to set the origins card of the current object
              * Set the origins_ attribute of the current pdb file
              * @param origins The origins attribute of the current object
              */
            void SetOrigins(PdbOriginXnSection* origins);
            /*! \fn
              * A mutator function in order to set the scales card of the current object
              * Set the scales_ attribute of the current pdb file
              * @param scales The scales attribute of the current object
              */
            void SetScales(PdbScaleNSection* scales);
            /*! \fn
              * A mutator function in order to set the matrices card of the current object
              * Set the matrices_ attribute of the current pdb file
              * @param matrices The matrices attribute of the current object
              */
            void SetMatrices(PdbMatrixNSection* matrices);
            /*! \fn
              * A mutator function in order to set the models card of the current object
              * Set the models_ attribute of the current pdb file
              * @param models The models attribute of the current object
              */
            void SetModels(PdbModelSection* models);
            /*! \fn
              * A mutator function in order to set the connectivities card of the current object
              * Set the connectivities_ attribute of the current pdb file
              * @param connectivities The connectivities attribute of the current object
              */
            void SetConnectivities(PdbConnectSection* connectivities);
            /*! \fn
              * A mutator function in order to set the serial number mapping of the current object
              * Set the serial_number_mapping_ attribute of the current pdb file
              * @param serial_number_mapping The serial number mapping attribute of the current object
              */
            void SetSerialNumberMapping(PdbSerialNumberMapping serial_number_mapping);
            /*! \fn
              * A mutator function in order to set the sequence number mapping of the current object
              * Set the sequence_number_mapping_ attribute of the current pdb file
              * @param sequence_number_mapping The sequence number mapping attribute of the current object
              */
            void SetSequenceNumberMapping(PdbSequenceNumberMapping sequence_number_mapping);
            /*! \fn
              * A function in order to delete a residue from the current object
              * @param residue A residue of the current object of this class
              */
            void DeleteResidue(PdbResidue* residue);
            /*! \fn
              * A function in order to delete a list of residue from the current object
              * @param residues List of residues of the current object of this class
              */
            void DeleteResidues(PdbResidueVector residues);
            /*! \fn
              * A function in order to delete a residue from the current object
              * @param residue A residue of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void DeleteResidueWithTheGivenModelNumber(PdbResidue* residue, int model_number = 1);
            /*! \fn
              * A function in order to delete a list of residues from the current object
              * @param residues List of residues of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void DeleteResiduesWithTheGivenModelNumber(PdbResidueVector residues, int model_number = 1);
            /*! \fn
              * A function in order to delete an atom from the current object
              * @param atom An atom of the current object of this class
              */
            void DeleteAtom(PdbAtomCard* atom);
            /*! \fn
              * A function in order to delete a list of atoms from the current object
              * @param atoms List of atoms of the current object of this class
              */
            void DeleteAtoms(PdbAtomCardVector atoms);
            /*! \fn
              * A function in order to delete an atom from the current object
              * @param atom An atom of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void DeleteAtomWithTheGivenModelNumber(PdbAtomCard* atom, int model_number = 1);
            /*! \fn
              * A function in order to delete a list atoms from the current object
              * @param atoms List of atoms of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void DeleteAtomsWithTheGivenModelNumber(PdbAtomCardVector atoms, int model_number = 1);
            /*! \fn
              * A function in order to update the residue name of a residue of the current object
              * @param residue A residue of the current object of this class
              * @param residue_name The new residue name for the given residue of the current object of this class
              */
            void UpdateResidueName(PdbResidue* residue, std::string residue_name);
            /*! \fn
              * A function in order to update the residue name of a residue of the current object
              * @param residue A residue of the current object of this class
              * @param residue_name The new residue name for the given residue of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void UpdateResidueNameWithTheGivenModelNumber(PdbResidue* residue, std::string residue_name, int model_number = 1);
            /*! \fn
              * A function in order to insert a residue before the given residue in a chain
              * @param residue A residue of the current object of this class
              */
            void InsertResidueBefore(PdbAtomSection* residue);
            /*! \fn
              * A function in order to insert a residue before the given residue in a chain
              * @param residue A residue of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void InsertResidueBeforeWithTheGivenModelNumber(PdbAtomSection* residue, int model_number = 1);
            /*! \fn
              * A function in order to insert a residue after the given residue in a chain
              * @param residue A residue of the current object of this class
              */
            void InsertResidueAfter(PdbAtomSection* residue);
            /*! \fn
              * A function in order to insert a residue after the given residue in a chain
              * @param residue A residue of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void InsertResidueAfterWithTheGivenModelNumber(PdbAtomSection* residue, int model_number = 1);
            /*! \fn
              * A function in order to split a chian of a model in atom card and put ter card at the given point
              * @param split_point_chain_id Residue chain id at the split point
              * @param split_point_sequence_number
              */
            void SplitAtomCardOfModelCard(char split_point_chain_id, int split_point_sequence_number);
            /*! \fn
              * A function in order to split a chian of a model in atom card and put ter card at the given point
              * @param split_point_chain_id Residue chain id at the split point
              * @param split_point_sequence_number
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void SplitAtomCardOfModelCardWithTheGivenModelNumber(char split_point_chain_id, int split_point_sequence_number, int model_number = 1);
            /*! \fn
              * A function to adjust the serial numbers of atoms that have been changed in connect card
              */
            void UpdateConnectCard();
            /*! \fn
              * An accessor function in order to set the master card of the current object of this class
              * @param master The master card of the current object
              * @return master_ The master card of the current object of this class
              */
            void SetMasterCard(PdbMasterCard* master);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Input_File_Reader
              * @{
              */
            /// Readers

            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a pdb file
              */
            bool Read(std::ifstream& in_file);
            /*! \fn
              * A function to parse the contents of a given stringstream of atoms
              * Parse the given stream and set the attributes of the current object accordingly
              * @param atomstream A stringstream contains atoms cards of a pdb file
              */
            bool Read(std::stringstream& atomstream);
            /*! \fn
              * A function to parse the contents of a given stringstream of atoms
              * @param atomstream A stream contains atom cards of a pdb file
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseAtomStream(std::stringstream& atomstream);
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * @param in_stream A stream contains whole contents of a pdb file
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCards(std::ifstream& in_stream);
            /*! \fn
              * A function to parse the header card that has been given as a stream
              * @param stream A stream contains header card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHeaderCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the obsolete card that has been given as a stream
              * @param stream A stream contains obsolete card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseObsoleteSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the title card that has been given as a stream
              * @param stream A stream contains title card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseTitleSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the split card that has been given as a stream
              * @param stream A stream contains split card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSplitSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the caveat card that has been given as a stream
              * @param stream A stream contains caveat card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCaveatSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the compound card that has been given as a stream
              * @param stream A stream contains compound card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCompoundSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the source card that has been given as a stream
              * @param stream A stream contains source card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSourceSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the keyword card that has been given as a stream
              * @param stream A stream contains keyword card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseKeywordsSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the expiration date card that has been given as a stream
              * @param stream A stream contains expiration date card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseExperimentalDataSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the num model card that has been given as a stream
              * @param stream A stream contains num model card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseNumModelCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the model type card that has been given as a stream
              * @param stream A stream contains model type card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseModelTypeSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the author card that has been given as a stream
              * @param stream A stream contains author card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseAuthorSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the revision date card that has been given as a stream
              * @param stream A stream contains revision date card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseRevisionDataSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the superseded entries card that has been given as a stream
              * @param stream A stream contains superseded entries card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSupersededEntriesSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the journal card that has been given as a stream
              * @param stream A stream contains journal card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseJournalSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the remark card that has been given as a stream
              * @param stream A stream contains remark card of a pdb file
              * @param lHeaderine Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseRemarkSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the database reference card that has been given as a stream
              * @param stream A stream contains database reference card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseDatabaseReferenceSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the sequence advanced card that has been given as a stream
              * @param stream A stream contains sequence advanced card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSequenceAdvancedSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the sequence residue card that has been given as a stream
              * @param stream A stream contains sequence reisdue card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseResidueSequenceSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the modification residue card that has been given as a stream
              * @param stream A stream contains modification residue card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseResidueModificationSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the heterogen card that has been given as a stream
              * @param stream A stream contains heterogen card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHeterogenSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the heterogen name card that has been given as a stream
              * @param stream A stream contains heterogen name card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHeterogenNameSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the heterogen synonym card that has been given as a stream
              * @param stream A stream contains heterogen synonym card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHeterogenSynonymSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the formula card that has been given as a stream
              * @param stream A stream contains formula card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseFormulaSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the helix card that has been given as a stream
              * @param stream A stream contains helix card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHelixSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the sheet card that has been given as a stream
              * @param stream A stream contains sheet card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSheetSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the disulfide bond card that has been given as a stream
              * @param stream A stream contains disulfide bond card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseDisulfideBondSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the link card that has been given as a stream
              * @param stream A stream contains link card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseLinkSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the cis peptide card that has been given as a stream
              * @param stream A stream contains cis peptide card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCISPeptideSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the site card that has been given as a stream
              * @param stream A stream contains site card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSiteSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the crystallography card that has been given as a stream
              * @param stream A stream contains crystallography card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCrystallographyCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the origin card that has been given as a stream
              * @param stream A stream contains origin card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseOriginCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the scale card that has been given as a stream
              * @param stream A stream contains scale card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseScaleCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the matrix card that has been given as a stream
              * @param stream A stream contains matrix card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseMatrixSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the model card that has been given as a stream
              * @param stream A stream contains model card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseModelSection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the connectivity card that has been given as a stream
              * @param stream A stream contains connectivity card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseConnectivitySection(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the master card that has been given as a stream
              * @param stream A stream contains master card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseMasterCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the end card that has been given as a stream
              * @param stream A stream contains end card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseEndCard(std::ifstream& stream, std::string& line);
/** @}*/
            /** \addtogroup Output_File_Builder
              * @{
              */
            /// Writers
            /*! \fn
              * A function to create an output pdb file in a stringstream
              * @param  pdbstream Output pdb stringstream
              */
            void WriteToStringstream(std::ostringstream& pdbstream);
            /*! \fn
              * A function to create an output pdb file with the given name
              * @param pdb_file Output pdb file name
              */
            void Write(const std::string& pdb_file);
            /*! \fn
              * A function to create an output pdb file with the given name
              * @param pdb_file Output pdb file name
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void WriteWithTheGivenModelNumber(const std::string& pdb_file, int model_number = 1);
            /*! \fn
              * A function to write back all card of the pdb file into an output stream
              * @param out_stream Primary output stream to write into an output file
              */
            void ResolveCards(std::ostream& out_stream);
            /*! \fn
              * A function to write back all card of the pdb file into an output stream
              * @param out_stream Primary output stream to write into an output file
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void ResolveCardsWithTheGivenModelNumber(std::ofstream& out_stream, int model_number = 1);
            /*! \fn
              * A function to write back header card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write header card
              */
            void ResolveHeaderCard(std::ostream& stream);
            /*! \fn
              * A function to write back obsolete card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write obsolete card
              */
            void ResolveObsoleteCards(std::ostream& stream);
            /*! \fn
              * A function to write back title card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write title card
              */
            void ResolveTitleCards(std::ostream& stream);
            /*! \fn
              * A function to write back split card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write split card
              */
            void ResolveSplitCards(std::ostream& stream);
            /*! \fn
              * A function to write back caveat card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write caveat card
              */
            void ResolveCaveatCards(std::ostream& stream);
            /*! \fn
              * A function to write back compound card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write compound card
              */
            void ResolveCompoundCards(std::ostream& stream);
            /*! \fn
              * A function to write back source card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write source card
              */
            void ResolveSourceCards(std::ostream& stream);
            /*! \fn
              * A function to write back keyword card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write keyword card
              */
            void ResolveKeywordCards(std::ostream& stream);
            /*! \fn
              * A function to write back expiration date card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write expiration date card
              */
            void ResolveExperimentalDataCards(std::ostream& stream);
            /*! \fn
              * A function to write back num model card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write num model card
              */
            void ResolveNumModelCard(std::ostream& stream);
            /*! \fn
              * A function to write back model type card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write model type card
              */
            void ResolveModelTypeCards(std::ostream& stream);
            /*! \fn
              * A function to write back author card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write author card
              */
            void ResolveAuthorCards(std::ostream& stream);
            /*! \fn
              * A function to write back revision date card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write revision date card
              */
            void ResolveRevisionDataCards(std::ostream& stream);
            /*! \fn
              * A function to write back superseded entries card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write superseded entries card
              */
            void ResolveSupersededEntriesCards(std::ostream& stream);
            /*! \fn
              * A function to write back journal card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write journal card
              */
            void ResolveJournalCards(std::ostream& stream);
            /*! \fn
              * A function to write back remark card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write remark card
              */
            void ResolveRemarkCards(std::ostream& stream);
            /*! \fn
              * A function to write back database reference card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write database reference card
              */
            void ResolveDatabaseReferenceCards(std::ostream& stream);
            /*! \fn
              * A function to write back sequence advanced card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write sequence advanced card
              */
            void ResolveSequenceAdvancedCards(std::ostream& stream);
            /*! \fn
              * A function to write back sequence residue card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write sequence residue card
              */
            void ResolveSequenceResidueCards(std::ostream& stream);
            /*! \fn
              * A function to write back modification residue card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write modification residue card
              */
            void ResolveModificationResidueCards(std::ostream& stream);
            /*! \fn
              * A function to write back heterogen card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write heterogen card
              */
            void ResolveHeterogenCards(std::ostream& stream);
            /*! \fn
              * A function to write back heterogen name card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write heterogen name card
              */
            void ResolveHeterogenNameCards(std::ostream& stream);
            /*! \fn
              * A function to write back heterogen synonym card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write heterogen synonym card
              */
            void ResolveHeterogenSynonymCards(std::ostream& stream);
            /*! \fn
              * A function to write back formula card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write formula card
              */
            void ResolveFormulaCards(std::ostream& stream);
            /*! \fn
              * A function to write back helix card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write helix card
              */
            void ResolveHelixCards(std::ostream& stream);
            /*! \fn
              * A function to write back sheet card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write sheet card
              */
            void ResolveSheetCards(std::ostream& stream);
            /*! \fn
              * A function to write back disulfide bond card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write disulfide bond card
              */
            void ResolveDisulfideBondCards(std::ostream& stream);
            /*! \fn
              * A function to write back link card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write link card
              */
            void ResolveLinkCards(std::ostream& stream);
            /*! \fn
              * A function to write back cis peptide card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write cis peptide card
              */
            void ResolveCISPeptideCards(std::ostream& stream);
            /*! \fn
              * A function to write back site card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write site card
              */
            void ResolveSiteCards(std::ostream& stream);
            /*! \fn
              * A function to write back crystallography card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write crystallography card
              */
            void ResolveCrystallographyCard(std::ostream& stream);
            /*! \fn
              * A function to write back origin card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write origin card
              */
            void ResolveOriginCard(std::ostream& stream);
            /*! \fn
              * A function to write back scale card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write scale card
              */
            void ResolveScaleCard(std::ostream& stream);
            /*! \fn
              * A function to write back matrix card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write matrix card
              */
            void ResolveMatrixCards(std::ostream& stream);
            /*! \fn
              * A function to write back model card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write model card
              */
            void ResolveModelCards(std::ostream& stream);
            /*! \fn
              * A function to write back model card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write model card
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void ResolveModelCardWithTheGivenModelNumber(std::ostream& stream, int model_number = 1);
            /*! \fn
              * A function to write back connectivity card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write connectivity card
              */
            void ResolveConnectivityCards(std::ostream& stream);
            /*! \fn
              * A function to write back master card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write master card
              */
            void ResolveMasterCards(std::ostream& stream);
            /*! \fn
              * A function to write back end card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write end card
              */
            void ResolveEndCard(std::ostream& stream);
            /*! \fn
              * A function to write relevant bits of pdb file into an output stream formatted for the ontology
              * @param stream output stream in order to write ontology
              */
            void PrintOntology(std::stringstream& ont_stream);
/** @}*/
            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string path_;                                          /*!< Path of the given pdb file >*/
            PdbHeaderCard* header_;                                     /*!< Header card >*/
            PdbObsoleteSection* obsolete_;                              /*!<Obsolete section>*/
            PdbTitleSection* title_;                                    /*!< Title section >*/
            PdbSplitSection* split_;                                    /*!< Split Section >*/
            PdbCaveatSection* caveat_;                                  /*!< Caveat Section>*/
            PdbCompoundSection* compound_;                              /*!< Compound section >*/
            PdbSourceSection* source_;                                  /*!< Source section >*/
            PdbKeywordsSection* keywords_;                               /*!< Keywords section >*/
            PdbExperimentalDataSection* experimental_data_;             /*!< Experimental data section>*/
            PdbNumModelCard* number_of_models_;                         /*!< Number of models card >*/
            PdbModelTypeSection* model_type_;                           /*!< Model type section >*/
            PdbAuthorSection* author_;                                  /*!< Author section >*/
            PdbRevisionDataSection* revision_data_;                     /*!< Revision data section >*/
            PdbSupersededEntriesSection* superseded_entries_;           /*!< Superseded entries section >*/
            PdbJournalSection* journal_;                                /*!< Journal section >*/
            PdbRemarkSection* remark_cards_;                            /*!< Remarks section>*/
            PdbDatabaseReferenceSection* database_reference_;           /*!< Database reference section >*/
            PdbSequenceAdvancedSection* sequence_advanced_;             /*!< Sequence advanced section >*/
            PdbResidueSequenceSection* residues_sequence_;              /*!< Residue sequence section >*/
            PdbResidueModificationSection* residue_modification_cards_; /*!< Residue modification section >*/
            PdbHeterogenSection* heterogen_cards_;                      /*!< Heterogen section >*/
            PdbHeterogenNameSection* heterogen_name_cards_;             /*!< Heterogen name section >*/
            PdbHeterogenSynonymSection* heterogen_synonym_cards_;       /*!< Heterogen synonym section >*/
            PdbFormulaSection* formulas_;                               /*!< Formula section >*/
            PdbHelixSection* helix_cards_;                              /*!< Helix section >*/
            PdbSheetSection* sheet_cards_;                              /*!< Sheet section >*/
            PdbDisulfideBondSection* disulfide_bonds_;                  /*!< Disulfide bond section >*/
            PdbLinkSection* link_cards_;                                /*!< Link section >*/
            PdbCISPeptideSection* cis_peptide_;                         /*!< CIS peptide section >*/
            PdbSiteSection* site_cards_;                                /*!< Site section >*/
            PdbCrystallographicCard* crystallography_;                  /*!< Crystallography card >*/
            PdbOriginXnSection* origins_;                               /*!< Origin section >*/
            PdbScaleNSection* scales_;                                  /*!< Scale section >*/
            PdbMatrixNSection* matrices_;                               /*!< Matrix section >*/
            PdbModelSection* models_;                                   /*!< Model section >*/
            PdbConnectSection* connectivities_;                         /*!< Connectivity section >*/
            PdbSerialNumberMapping serial_number_mapping_;              /*!< A map that keeps track of serial numbers that have been changed during a process >*/
            PdbSequenceNumberMapping sequence_number_mapping_;          /*!< A map that keeps track of sequence numbers that have been changed during a process >*/
            PdbMasterCard* master_;                                     /*!< Master card >*/
            // PdbEndCard* end_;                                           /*!< End Card >*/
    };
}

#endif // PDBFILE_HPP
