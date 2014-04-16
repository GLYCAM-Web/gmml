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

    class PdbFile
    {
        public:            
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<PdbResidue*> PdbResidueVector;
            typedef std::vector<PdbAtom*> PdbAtomVector;
            typedef std::map<std::string, PdbAtomVector* > PdbResidueAtomsMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor
              * @param pdb_file An existing pdb file path to be read
              */
            PdbFile(const std::string& pdb_file);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetPath();
            PdbHeaderCard* GetHeader();
            PdbTitleCard* GetTitle();
            PdbCompoundCard* GetCompound();
            PdbNumModelCard* GetNumberOfModels();
            PdbModelTypeCard* GetModelType();
            PdbResidueSequenceCard* GetResiduesSequence();
            PdbResidueModificationCard* GetResidueModification();
            PdbHeterogenCard* GetHeterogens();
            PdbHeterogenNameCard* GetHeterogensName();
            PdbHeterogenSynonymCard* GetHeterogenSynonyms();
            PdbFormulaCard* GetFormulas();
            PdbHelixCard* GetHelixes();
            PdbSheetCard* GetSheets();
            PdbDisulfideBondCard* GetDisulfideBonds();
            PdbLinkCard* GetLinks();
            PdbSiteCard* GetSites();
            PdbCrystallographicCard* GetCrystallography();
            PdbOriginXnCard* GetOrigins();
            PdbScaleNCard* GetScales();
            PdbMatrixNCard* GetMatrices();
            PdbModelCard* GetModels();
            PdbConnectCard* GetConnectivities();            
            std::vector<std::string> GetAllResidueNames();
            PdbResidueVector GetAllResidues();
            PdbAtomVector GetAllAtomsOfResidue(PdbResidue* residue);
            PdbResidueAtomsMap GetAllAtomsOfResidues();
            PdbAtom* GetAtomOfResidueByName(PdbResidue* residue, std::string atom_name, PdbResidueAtomsMap residue_atom_map);
            PdbAtom* GetAtomOfResidueByName(PdbResidue* residue, std::string atom_name);
            std::vector<std::string> GetAllAtomNamesOfResidue(PdbResidue* residue, PdbResidueAtomsMap residue_atom_map);
            std::vector<std::string> GetAllAtomNamesOfResidue(PdbResidue* residue);

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////

            void DeleteResidue(PdbResidue* residue);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            /// Readers

            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a pdb file
              */
            void Read(std::ifstream& in_file);
            void ParseCards(std::ifstream& in_stream);
            void ParseHeaderCard(std::ifstream& stream, std::string& line);
            void ParseObsoleteCard(std::ifstream& stream, std::string& line);
            void ParseTitleCard(std::ifstream& stream, std::string& line);
            void ParseSplitCard(std::ifstream& stream, std::string& line);
            void ParseCaveatCard(std::ifstream& stream, std::string& line);
            void ParseCompoundCard(std::ifstream& stream, std::string& line);
            void ParseSourceCard(std::ifstream& stream, std::string& line);
            void ParseKeywordCard(std::ifstream& stream, std::string& line);
            void ParseExpirationDateCard(std::ifstream& stream, std::string& line);
            void ParseNumModelCard(std::ifstream& stream, std::string& line);
            void ParseModelTypeCard(std::ifstream& stream, std::string& line);
            void ParseAuthorCard(std::ifstream& stream, std::string& line);
            void ParseRevisionDateCard(std::ifstream& stream, std::string& line);
            void ParseSupersededEntriesCard(std::ifstream& stream, std::string& line);
            void ParseJournalCard(std::ifstream& stream, std::string& line);
            void ParseRemarkCard(std::ifstream& stream, std::string& line);
            void ParseDatabaseReferenceCard(std::ifstream& stream, std::string& line);
            void ParseSequenceAdvancedCard(std::ifstream& stream, std::string& line);
            void ParseSequenceResidueCard(std::ifstream& stream, std::string& line);
            void ParseModificationResidueCard(std::ifstream& stream, std::string& line);
            void ParseHeterogenCard(std::ifstream& stream, std::string& line);
            void ParseHeterogenNameCard(std::ifstream& stream, std::string& line);
            void ParseHeterogenSynonymCard(std::ifstream& stream, std::string& line);
            void ParseFormulaCard(std::ifstream& stream, std::string& line);
            void ParseHelixCard(std::ifstream& stream, std::string& line);
            void ParseSheetCard(std::ifstream& stream, std::string& line);
            void ParseDisulfideBondCard(std::ifstream& stream, std::string& line);
            void ParseLinkCard(std::ifstream& stream, std::string& line);
            void ParseCISPeptideCard(std::ifstream& stream, std::string& line);
            void ParseSiteCard(std::ifstream& stream, std::string& line);
            void ParseCrystallographyCard(std::ifstream& stream, std::string& line);
            void ParseOriginCard(std::ifstream& stream, std::string& line);
            void ParseScaleCard(std::ifstream& stream, std::string& line);
            void ParseMatrixCard(std::ifstream& stream, std::string& line);
            void ParseModelCard(std::ifstream& stream, std::string& line);
            void ParseConnectivityCard(std::ifstream& stream, std::string& line);
            void ParseMasterCard(std::ifstream& stream, std::string& line);
            void ParseEndCard(std::ifstream& stream, std::string& line);

            /// Writers
            void Write(const std::string& pdb_file);
            void ResolveCards(std::ofstream& out_stream);
            void ResolveHeaderCard(std::ofstream& stream);
            void ResolveObsoleteCard(std::ofstream& stream);
            void ResolveTitleCard(std::ofstream& stream);
            void ResolveSplitCard(std::ofstream& stream);
            void CaveatCard(std::ofstream& stream);
            void ResolveCompoundCard(std::ofstream& stream);
            void ResolveSourceCard(std::ofstream& stream);
            void ResolveKeywordCard(std::ofstream& stream);
            void ResolveExpirationDateCard(std::ofstream& stream);
            void ResolveNumModelCard(std::ofstream& stream);
            void ResolveModelTypeCard(std::ofstream& stream);
            void ResolveAuthorCard(std::ofstream& stream);
            void ResolveRevisionDateCard(std::ofstream& stream);
            void ResolveSupersededEntriesCard(std::ofstream& stream);
            void ResolveJournalCard(std::ofstream& stream);
            void ResolveRemarkCard(std::ofstream& stream);
            void ResolveDatabaseReferenceCard(std::ofstream& stream);
            void ResolveSequenceAdvancedCard(std::ofstream& stream);
            void ResolveSequenceResidueCard(std::ofstream& stream);
            void ResolveModificationResidueCard(std::ofstream& stream);
            void ResolveHeterogenCard(std::ofstream& stream);
            void ResolveHeterogenNameCard(std::ofstream& stream);
            void ResolveHeterogenSynonymCard(std::ofstream& stream);
            void ResolveFormulaCard(std::ofstream& stream);
            void ResolveHelixCard(std::ofstream& stream);
            void ResolveSheetCard(std::ofstream& stream);
            void ResolveDisulfideBondCard(std::ofstream& stream);
            void ResolveLinkCard(std::ofstream& stream);
            void ResolveCISPeptideCard(std::ofstream& stream);
            void ResolveSiteCard(std::ofstream& stream);
            void ResolveCrystallographyCard(std::ofstream& stream);
            void ResolveOriginCard(std::ofstream& stream);
            void ResolveScaleCard(std::ofstream& stream);
            void ResolveMatrixCard(std::ofstream& stream);
            void ResolveModelCard(std::ofstream& stream);
            void ResolveConnectivityCard(std::ofstream& stream);
            void ResolveMasterCard(std::ofstream& stream);
            void ResolveEndCard(std::ofstream& stream);


            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string path_;
            PdbHeaderCard* header_;
            PdbTitleCard* title_;
            PdbCompoundCard* compound_;
            PdbNumModelCard* number_of_models_;
            PdbModelTypeCard* model_type_;
            PdbResidueSequenceCard* residues_sequence_;
            PdbResidueModificationCard* residue_modification_;
            PdbHeterogenCard* heterogens_;
            PdbHeterogenNameCard* heterogens_name_;
            PdbHeterogenSynonymCard* heterogen_synonyms_;
            PdbFormulaCard* formulas_;
            PdbHelixCard* helixes_;
            PdbSheetCard* sheets_;
            PdbDisulfideBondCard* disulfide_bonds_;
            PdbLinkCard* links_;
            PdbSiteCard* sites_;
            PdbCrystallographicCard* crystallography_;
            PdbOriginXnCard* origins_;
            PdbScaleNCard* scales_;
            PdbMatrixNCard* matrices_;
            PdbModelCard* models_;
            PdbConnectCard* connectivities_;
    };
}

#endif // PDBFILE_HPP
