//Author: Alireza Khatamian

#ifndef PDBFILE_HPP
#define PDBFILE_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace PdbFileSpace
{
    class PdbHeaderCard;
    class PdbTitleCard;
    class PdbCompoundCard;
    class PdbNumModelCard;
    class PdbModelTypeCard;
    class PdbResidueSequenceCard;
    class PdbResidueModificationCard;
    class PdbHeterogenCard;
    class PdbHeterogenNameCard;
    class PdbHeterogenSynonymCard;
    class PdbFormulaCard;
    class PdbHelixCard;
    class PdbSheetCard;
    class PdbDisulfideBondCard;
    class PdbLinkCard;
    class PdbSiteCard;
    class PdbCrystallographicCard;
    class PdbOriginXnCard;
    class PdbScaleNCard;
    class PdbMatrixNCard;
    class PdbModelCard;
    class PdbConnectCard;
    class PdbResidue;
    class PdbAtom;
    class PdbAtomCard;
    class PdbFile
    {
        public:            
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of residues
              */
            typedef std::vector<PdbResidue*> PdbResidueVector;
            /*! \typedef
              * List of pdb atom
              */
            typedef std::vector<PdbAtom*> PdbAtomVector;
            /*! \typedef
              * A mapping between a
              */
            typedef std::map<std::string, PdbAtomVector* > PdbResidueAtomsMap;
            /*! \typedef
              * Mapping between old serial number and new one that has been changed during a process
              */
            typedef std::map<int, int> PdbSerialNumberMapping;

            typedef std::vector<std::pair<char, int> > PdbPairVectorTerCardPositions;

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

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
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
              * An accessor function in order to access to the title of the current object
              * @return title_ attribute of the current object of this class
              */
            PdbTitleCard* GetTitle();
            /*! \fn
              * An accessor function in order to access to the compound attribute of the current object
              * @return compound_ attribute of the current object of this class
              */
            PdbCompoundCard* GetCompound();
            /*! \fn
              * An accessor function in order to access to the number of models attribute of the current object
              * @return number_of_models_ attribute of the current object of this class
              */
            PdbNumModelCard* GetNumberOfModels();
            /*! \fn
              * An accessor function in order to access to the model type attribute of the current object
              * @return model_type_ attribute of the current object of this class
              */
            PdbModelTypeCard* GetModelType();
            /*! \fn
              * An accessor function in order to access to the residue sequence attribute of the current object
              * @return residue_sequence_ attribute of the current object of this class
              */
            PdbResidueSequenceCard* GetResiduesSequence();
            /*! \fn
              * An accessor function in order to access to the residue modification attribute of the current object
              * @return residue_modification_ attribute of the current object of this class
              */
            PdbResidueModificationCard* GetResidueModification();
            /*! \fn
              * An accessor function in order to access to the heterogens attribute of the current object
              * @return heterogens_ attribute of the current object of this class
              */
            PdbHeterogenCard* GetHeterogens();
            /*! \fn
              * An accessor function in order to access to the heterogens name attribute of the current object
              * @return heterogens_name_ attribute of the current object of this class
              */
            PdbHeterogenNameCard* GetHeterogensName();
            /*! \fn
              * An accessor function in order to access to the heterogen synonyms attribute of the current object
              * @return heterogen_synonyms_ attribute of the current object of this class
              */
            PdbHeterogenSynonymCard* GetHeterogenSynonyms();
            /*! \fn
              * An accessor function in order to access to the formulas attribute of the current object
              * @return formulas_ attribute of the current object of this class
              */
            PdbFormulaCard* GetFormulas();
            /*! \fn
              * An accessor function in order to access to the helixes attribute of the current object
              * @return helixes_ attribute of the current object of this class
              */
            PdbHelixCard* GetHelixes();
            /*! \fn
              * An accessor function in order to access to the sheets attribute of the current object
              * @return sheets_ attribute of the current object of this class
              */
            PdbSheetCard* GetSheets();
            /*! \fn
              * An accessor function in order to access to the disulfide bonds attribute of the current object
              * @return disulfide_bonds_ attribute of the current object of this class
              */
            PdbDisulfideBondCard* GetDisulfideBonds();
            /*! \fn
              * An accessor function in order to access to the links attribute of the current object
              * @return links_ attribute of the current object of this class
              */
            PdbLinkCard* GetLinks();
            /*! \fn
              * An accessor function in order to access to the sites attribute of the current object
              * @return sites_ attribute of the current object of this class
              */
            PdbSiteCard* GetSites();
            /*! \fn
              * An accessor function in order to access to the crystallography attribute of the current object
              * @return crystallography_ attribute of the current object of this class
              */
            PdbCrystallographicCard* GetCrystallography();
            /*! \fn
              * An accessor function in order to access to the origins attribute of the current object
              * @return origins_ attribute of the current object of this class
              */
            PdbOriginXnCard* GetOrigins();
            /*! \fn
              * An accessor function in order to access to the scales attribute of the current object
              * @return scales_ attribute of the current object of this class
              */
            PdbScaleNCard* GetScales();
            /*! \fn
              * An accessor function in order to access to the matrices attribute of the current object
              * @return matrices_ attribute of the current object of this class
              */
            PdbMatrixNCard* GetMatrices();
            /*! \fn
              * An accessor function in order to access to the models attribute of the current object
              * @return models_ attribute of the current object of this class
              */
            PdbModelCard* GetModels();
            /*! \fn
              * An accessor function in order to access to the connectivities attribute of the current object
              * @return connectivities_ attribute of the current object of this class
              */
            PdbConnectCard* GetConnectivities();
            /*! \fn
              * An accessor function in order to access to the serial number mapping attribute of the current object
              * @return serial_number_mapping attribute of the current object of this class
              */
            PdbSerialNumberMapping GetSerialNumberMapping();
            /*! \fn
              * An accessor function in order to access to all residue names of the current object
              * @return residue_names All residue names of the current object of this class
              */
            std::vector<std::string> GetAllResidueNames();
            /*! \fn
              * An accessor function in order to access to all residue names from atom card of the current object
              * @return residue_names All residue names from atom card of the current object of this class
              */
            std::vector<std::string> GetAllResidueNamesFromAtomCard();
            /*! \fn
              * An accessor function in order to access to all residues of the current object
              * @return all_residues_ All residues of the current object of this class
              */
            PdbResidueVector GetAllResidues();
            /*! \fn
              * An accessor function in order to access to all residues from atom card of the current object
              * @return residues All resdidues from atom card of the current object of this class
              */
            PdbResidueVector GetAllResiduesFromAtomCard();
            /*! \fn
              * An accessor function in order to access to all atoms of a residue of the current object
              * @param residue The given residue to return all the atoms of it
              * @return atoms_of_residue All atoms of a resdidue of the current object of this class
              */
            PdbAtomVector GetAllAtomsOfResidue(PdbResidue* residue);
            /*! \fn
              * An accessor function in order to access to all atoms of all residues of the current object
              * @return residue_atom_map The map between all residues and their atoms of the current object of this class
              */
            PdbResidueAtomsMap GetAllAtomsOfResidues();
            /*! \fn
              * An accessor function in order to access to atom of a residue of the current object using the atom name
              * @param residue The given residue to return one of its atoms
              * @param atom_name The atom name of the desired atom object
              * @param residue_atom_map The map between residues and their atoms
              * @return atom The atom object of the given residue of the current object of this class
              */
            PdbAtom* GetAtomOfResidueByName(PdbResidue* residue, std::string atom_name, PdbResidueAtomsMap residue_atom_map);
            /*! \fn
              * An accessor function in order to access to atom of a residue of the current object using the atom name
              * @param residue The given residue to return one of its atoms
              * @param atom_name The atom name of the desired atom object
              * @return atom The atom object of the given residue of the current object of this class
              */
            PdbAtom* GetAtomOfResidueByName(PdbResidue* residue, std::string atom_name);
            /*! \fn
              * An accessor function in order to access to atom of a residue of the current object using the atom key
              * @param atom_key The atom key of the desired atom object
              * @return heterogen_atom_map The map between heterogen atoms of the current object of this class and their keys
              */
            PdbAtom* GetAtomOfResidueByAtomKey(std::string atom_key);
            /*! \fn
              * An accessor function in order to access to atom of a residue of the current object using the serial number
              * @param serial_number The serial number attribute of the desired atom object
              * @return heterogen_atom_map The map between heterogen atoms of the current object of this class and their keys
              */
            PdbAtom* GetAtomBySerialNumber(int serial_number);
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

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
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
              * A mutator function in order to set the title card of the current object
              * Set the title_ attribute of the current pdb file
              * @param title The title attribute of the current object
              */
            void SetTitle(PdbTitleCard* title);
            /*! \fn
              * A mutator function in order to set the compound card of the current object
              * Set the compound_ attribute of the current pdb file
              * @param compound The compound attribute of the current object
              */
            void SetCompound(PdbCompoundCard* compound);
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
            void SetModelType(PdbModelTypeCard* model_type);
            /*! \fn
              * A mutator function in order to set the residues sequence card of the current object
              * Set the residues_sequence_ attribute of the current pdb file
              * @param residues_sequence The residues sequence attribute of the current object
              */
            void SetResiduesSequence(PdbResidueSequenceCard* residues_sequence);
            /*! \fn
              * A mutator function in order to set the residue modification card of the current object
              * Set the residue_modification_ attribute of the current pdb file
              * @param residue_modification The residue modification attribute of the current object
              */
            void SetResidueModification(PdbResidueModificationCard* residue_modification);
            /*! \fn
              * A mutator function in order to set the heterogens card of the current object
              * Set the heterogens_ attribute of the current pdb file
              * @param heterogens The heterogens attribute of the current object
              */
            void SetHeterogens(PdbHeterogenCard* heterogens);
            /*! \fn
              * A mutator function in order to set the heterogens name card of the current object
              * Set the heterogens_name_ attribute of the current pdb file
              * @param heterogens_name The heterogens name attribute of the current object
              */
            void SetHeterogensName(PdbHeterogenNameCard* heterogens_name);
            /*! \fn
              * A mutator function in order to set the heterogen synonyms card of the current object
              * Set the heterogen_synonyms_ attribute of the current pdb file
              * @param heterogen_synonyms The heterogen synonyms attribute of the current object
              */
            void SetHeterogenSynonyms(PdbHeterogenSynonymCard* heterogen_synonyms);
            /*! \fn
              * A mutator function in order to set the formulas card of the current object
              * Set the formulas_ attribute of the current pdb file
              * @param formulas The formulas attribute of the current object
              */
            void SetFormulas(PdbFormulaCard* formulas);
            /*! \fn
              * A mutator function in order to set the helixes card of the current object
              * Set the helixes_ attribute of the current pdb file
              * @param helixes The helixes attribute of the current object
              */
            void SetHelixes(PdbHelixCard* helixes);
            /*! \fn
              * A mutator function in order to set the sheets card of the current object
              * Set the sheets_ attribute of the current pdb file
              * @param sheets The sheets attribute of the current object
              */
            void SetSheets(PdbSheetCard* sheets);
            /*! \fn
              * A mutator function in order to set the disulfide bonds card of the current object
              * Set the disulfide_bonds_ attribute of the current pdb file
              * @param disulfide_bonds The disulfide bonds attribute of the current object
              */
            void SetDisulfideBonds(PdbDisulfideBondCard* disulfide_bonds);
            /*! \fn
              * A mutator function in order to set the links card of the current object
              * Set the links_ attribute of the current pdb file
              * @param links The links attribute of the current object
              */
            void SetLinks(PdbLinkCard* links);
            /*! \fn
              * A mutator function in order to set the sites card of the current object
              * Set the sites_ attribute of the current pdb file
              * @param sites The compound attribute of the current object
              */
            void SetSites(PdbSiteCard* sites);
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
            void SetOrigins(PdbOriginXnCard* origins);
            /*! \fn
              * A mutator function in order to set the scales card of the current object
              * Set the scales_ attribute of the current pdb file
              * @param scales The scales attribute of the current object
              */
            void SetScales(PdbScaleNCard* scales);
            /*! \fn
              * A mutator function in order to set the matrices card of the current object
              * Set the matrices_ attribute of the current pdb file
              * @param matrices The matrices attribute of the current object
              */
            void SetMatrices(PdbMatrixNCard* matrices);
            /*! \fn
              * A mutator function in order to set the models card of the current object
              * Set the models_ attribute of the current pdb file
              * @param models The models attribute of the current object
              */
            void SetModels(PdbModelCard* models);
            /*! \fn
              * A mutator function in order to set the connectivities card of the current object
              * Set the connectivities_ attribute of the current pdb file
              * @param connectivities The connectivities attribute of the current object
              */
            void SetConnectivities(PdbConnectCard* connectivities);
            /*! \fn
              * A mutator function in order to set the serial number mapping of the current object
              * Set the serial_number_mapping_ attribute of the current pdb file
              * @param serial_number_mapping The serial number mapping attribute of the current object
              */
            void SetSerialNumberMapping(PdbSerialNumberMapping serial_number_mapping);
            /*! \fn
              * A function in order to delete a residue from the current object
              * @param residue A residue of the current object of this class
              */
            void DeleteResidue(PdbResidue* residue);
            /*! \fn
              * A function in order to delete a residue from the current object
              * @param residue A residue of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void DeleteResidueWithTheGivenModelNumber(PdbResidue* residue, int model_number = 1);
            /*! \fn
              * A function in order to delete an atom from the current object
              * @param atom An atom of the current object of this class
              */
            void DeleteAtom(PdbAtom* atom);
            /*! \fn
              * A function in order to delete an atom from the current object
              * @param atom An atom of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void DeleteAtomWithTheGivenModelNumber(PdbAtom* atom, int model_number = 1);
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
            void InsertResidueBefore(PdbAtomCard* residue);
            /*! \fn
              * A function in order to insert a residue before the given residue in a chain
              * @param residue A residue of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void InsertResidueBeforeWithTheGivenModelNumber(PdbAtomCard* residue, int model_number = 1);
            /*! \fn
              * A function in order to insert a residue after the given residue in a chain
              * @param residue A residue of the current object of this class
              */
            void InsertResidueAfter(PdbAtomCard* residue);
            /*! \fn
              * A function in order to insert a residue after the given residue in a chain
              * @param residue A residue of the current object of this class
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void InsertResidueAfterWithTheGivenModelNumber(PdbAtomCard* residue, int model_number = 1);
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

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            /// Readers

            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a pdb file
              */
            bool Read(std::ifstream& in_file);
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * @param in_stream A stream contains whole contents of a pdb file
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCards(std::ifstream& in_stream);
            /*! \fn
              * A function to parse the header crad that has been given as a stream
              * @param stream A stream contains header card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHeaderCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the obsolete crad that has been given as a stream
              * @param stream A stream contains obsolete card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseObsoleteCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the title crad that has been given as a stream
              * @param stream A stream contains title card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseTitleCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the split crad that has been given as a stream
              * @param stream A stream contains split card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSplitCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the caveat crad that has been given as a stream
              * @param stream A stream contains caveat card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCaveatCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the compound crad that has been given as a stream
              * @param stream A stream contains compound card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCompoundCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the source crad that has been given as a stream
              * @param stream A stream contains source card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSourceCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the keyword crad that has been given as a stream
              * @param stream A stream contains keyword card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseKeywordCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the expiration date crad that has been given as a stream
              * @param stream A stream contains expiration date card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseExpirationDateCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the num model crad that has been given as a stream
              * @param stream A stream contains num model card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseNumModelCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the model type crad that has been given as a stream
              * @param stream A stream contains model type card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseModelTypeCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the author crad that has been given as a stream
              * @param stream A stream contains author card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseAuthorCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the revision date crad that has been given as a stream
              * @param stream A stream contains revision date card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseRevisionDateCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the superseded entries crad that has been given as a stream
              * @param stream A stream contains superseded entries card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSupersededEntriesCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the journal crad that has been given as a stream
              * @param stream A stream contains journal card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseJournalCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the remark crad that has been given as a stream
              * @param stream A stream contains remark card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseRemarkCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the database reference crad that has been given as a stream
              * @param stream A stream contains database reference card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseDatabaseReferenceCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the sequence advanced crad that has been given as a stream
              * @param stream A stream contains sequence advanced card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSequenceAdvancedCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the sequence residue crad that has been given as a stream
              * @param stream A stream contains sequence reisdue card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSequenceResidueCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the modification residue crad that has been given as a stream
              * @param stream A stream contains modification residue card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseModificationResidueCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the heterogen crad that has been given as a stream
              * @param stream A stream contains heterogen card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHeterogenCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the heterogen name crad that has been given as a stream
              * @param stream A stream contains heterogen name card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHeterogenNameCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the heterogen synonym crad that has been given as a stream
              * @param stream A stream contains heterogen synonym card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHeterogenSynonymCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the formula crad that has been given as a stream
              * @param stream A stream contains formula card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseFormulaCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the helix crad that has been given as a stream
              * @param stream A stream contains helix card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseHelixCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the sheet crad that has been given as a stream
              * @param stream A stream contains sheet card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSheetCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the disulfide bond crad that has been given as a stream
              * @param stream A stream contains disulfide bond card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseDisulfideBondCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the link crad that has been given as a stream
              * @param stream A stream contains link card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseLinkCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the cis peptide crad that has been given as a stream
              * @param stream A stream contains cis peptide card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCISPeptideCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the site crad that has been given as a stream
              * @param stream A stream contains site card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseSiteCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the crystallography crad that has been given as a stream
              * @param stream A stream contains crystallography card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCrystallographyCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the origin crad that has been given as a stream
              * @param stream A stream contains origin card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseOriginCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the scale crad that has been given as a stream
              * @param stream A stream contains scale card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseScaleCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the matrix crad that has been given as a stream
              * @param stream A stream contains matrix card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseMatrixCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the model crad that has been given as a stream
              * @param stream A stream contains model card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseModelCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the connectivity crad that has been given as a stream
              * @param stream A stream contains connectivity card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseConnectivityCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the master crad that has been given as a stream
              * @param stream A stream contains master card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseMasterCard(std::ifstream& stream, std::string& line);
            /*! \fn
              * A function to parse the end crad that has been given as a stream
              * @param stream A stream contains end card of a pdb file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseEndCard(std::ifstream& stream, std::string& line);

            /// Writers
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
            void ResolveCards(std::ofstream& out_stream);
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
            void ResolveHeaderCard(std::ofstream& stream);
            /*! \fn
              * A function to write back obsolete card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write obsolete card
              */
            void ResolveObsoleteCard(std::ofstream& stream);
            /*! \fn
              * A function to write back title card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write title card
              */
            void ResolveTitleCard(std::ofstream& stream);
            /*! \fn
              * A function to write back split card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write split card
              */
            void ResolveSplitCard(std::ofstream& stream);
            /*! \fn
              * A function to write back caveat card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write caveat card
              */
            void ResolveCaveatCard(std::ofstream& stream);
            /*! \fn
              * A function to write back compound card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write compound card
              */
            void ResolveCompoundCard(std::ofstream& stream);
            /*! \fn
              * A function to write back source card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write source card
              */
            void ResolveSourceCard(std::ofstream& stream);
            /*! \fn
              * A function to write back keyword card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write keyword card
              */
            void ResolveKeywordCard(std::ofstream& stream);
            /*! \fn
              * A function to write back expiration date card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write expiration date card
              */
            void ResolveExpirationDateCard(std::ofstream& stream);
            /*! \fn
              * A function to write back num model card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write num model card
              */
            void ResolveNumModelCard(std::ofstream& stream);
            /*! \fn
              * A function to write back model type card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write model type card
              */
            void ResolveModelTypeCard(std::ofstream& stream);
            /*! \fn
              * A function to write back author card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write author card
              */
            void ResolveAuthorCard(std::ofstream& stream);
            /*! \fn
              * A function to write back revision date card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write revision date card
              */
            void ResolveRevisionDateCard(std::ofstream& stream);
            /*! \fn
              * A function to write back superseded entries card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write superseded entries card
              */
            void ResolveSupersededEntriesCard(std::ofstream& stream);
            /*! \fn
              * A function to write back journal card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write journal card
              */
            void ResolveJournalCard(std::ofstream& stream);
            /*! \fn
              * A function to write back remark card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write remark card
              */
            void ResolveRemarkCard(std::ofstream& stream);
            /*! \fn
              * A function to write back database reference card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write database reference card
              */
            void ResolveDatabaseReferenceCard(std::ofstream& stream);
            /*! \fn
              * A function to write back sequence advanced card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write sequence advanced card
              */
            void ResolveSequenceAdvancedCard(std::ofstream& stream);
            /*! \fn
              * A function to write back sequence residue card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write sequence residue card
              */
            void ResolveSequenceResidueCard(std::ofstream& stream);
            /*! \fn
              * A function to write back modification residue card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write modification residue card
              */
            void ResolveModificationResidueCard(std::ofstream& stream);
            /*! \fn
              * A function to write back heterogen card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write heterogen card
              */
            void ResolveHeterogenCard(std::ofstream& stream);
            /*! \fn
              * A function to write back heterogen name card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write heterogen name card
              */
            void ResolveHeterogenNameCard(std::ofstream& stream);
            /*! \fn
              * A function to write back heterogen synonym card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write heterogen synonym card
              */
            void ResolveHeterogenSynonymCard(std::ofstream& stream);
            /*! \fn
              * A function to write back formula card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write formula card
              */
            void ResolveFormulaCard(std::ofstream& stream);
            /*! \fn
              * A function to write back helix card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write helix card
              */
            void ResolveHelixCard(std::ofstream& stream);
            /*! \fn
              * A function to write back sheet card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write sheet card
              */
            void ResolveSheetCard(std::ofstream& stream);
            /*! \fn
              * A function to write back disulfide bond card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write disulfide bond card
              */
            void ResolveDisulfideBondCard(std::ofstream& stream);
            /*! \fn
              * A function to write back link card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write link card
              */
            void ResolveLinkCard(std::ofstream& stream);
            /*! \fn
              * A function to write back cis peptide card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write cis peptide card
              */
            void ResolveCISPeptideCard(std::ofstream& stream);
            /*! \fn
              * A function to write back site card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write site card
              */
            void ResolveSiteCard(std::ofstream& stream);
            /*! \fn
              * A function to write back crystallography card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write crystallography card
              */
            void ResolveCrystallographyCard(std::ofstream& stream);
            /*! \fn
              * A function to write back origin card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write origin card
              */
            void ResolveOriginCard(std::ofstream& stream);
            /*! \fn
              * A function to write back scale card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write scale card
              */
            void ResolveScaleCard(std::ofstream& stream);
            /*! \fn
              * A function to write back matrix card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write matrix card
              */
            void ResolveMatrixCard(std::ofstream& stream);
            /*! \fn
              * A function to write back model card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write model card
              */
            void ResolveModelCard(std::ofstream& stream);
            /*! \fn
              * A function to write back model card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write model card
              * @param model_number Selected model number from the multiple models that are in a pdb file
              */
            void ResolveModelCardWithTheGivenModelNumber(std::ofstream& stream, int model_number = 1);
            /*! \fn
              * A function to write back connectivity card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write connectivity card
              */
            void ResolveConnectivityCard(std::ofstream& stream);
            /*! \fn
              * A function to write back master card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write master card
              */
            void ResolveMasterCard(std::ofstream& stream);
            /*! \fn
              * A function to write back end card of the pdb file into an output stream
              * @param stream Intermediate output stream in order to write end card
              */
            void ResolveEndCard(std::ofstream& stream);

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string path_;                                      /*!< Path of the given pdb file >*/
            PdbHeaderCard* header_;                                 /*!< Header card >*/
            PdbTitleCard* title_;                                   /*!< Title carde >*/
            PdbCompoundCard* compound_;                             /*!< Compound card >*/
            PdbNumModelCard* number_of_models_;                     /*!< Number of models card >*/
            PdbModelTypeCard* model_type_;                          /*!< Model type card >*/
            PdbResidueSequenceCard* residues_sequence_;             /*!< Residue sequence card >*/
            PdbResidueModificationCard* residue_modification_;      /*!< Residue modification card >*/
            PdbHeterogenCard* heterogens_;                          /*!< Heterogen card >*/
            PdbHeterogenNameCard* heterogens_name_;                 /*!< Heterogen name card >*/
            PdbHeterogenSynonymCard* heterogen_synonyms_;           /*!< Heterogen synonym card >*/
            PdbFormulaCard* formulas_;                              /*!< Formula card >*/
            PdbHelixCard* helixes_;                                 /*!< Helix card >*/
            PdbSheetCard* sheets_;                                  /*!< Sheet card >*/
            PdbDisulfideBondCard* disulfide_bonds_;                 /*!< Disulfide bond card >*/
            PdbLinkCard* links_;                                    /*!< Link card >*/
            PdbSiteCard* sites_;                                    /*!< Site card >*/
            PdbCrystallographicCard* crystallography_;              /*!< Crystallography card >*/
            PdbOriginXnCard* origins_;                              /*!< Origin card >*/
            PdbScaleNCard* scales_;                                 /*!< Scale card >*/
            PdbMatrixNCard* matrices_;                              /*!< Matrix card >*/
            PdbModelCard* models_;                                  /*!< Model card >*/
            PdbConnectCard* connectivities_;                        /*!< Connectivity card >*/
            PdbSerialNumberMapping serial_number_mapping_;          /*!< A map that keeps track of serial numbers that have been changed during a process >*/
    };
}

#endif // PDBFILE_HPP
