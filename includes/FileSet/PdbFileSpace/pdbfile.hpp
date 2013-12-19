//Author: Alireza Khatamian

#ifndef PDBFILE_HPP
#define PDBFILE_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbHeaderCard;
    class PdbTitleCard;
    class PdbCompoundCard;
    class PdbNumModelCard;
    class PdbModelTypeCard;
    class PdbResidueSeqenceCard;
    class PdbResidueModificationCard;
    class PdbHeterogenCard;
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

    class PdbFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

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
            PdbResidueSeqenceCard* GetResidueSequence();
            PdbResidueModificationCard* GetResidueModification();
            PdbHeterogenCard* GetHeterogens();
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

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a pdb file
              */
            void Read(std::ifstream& in_file);

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////

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
            PdbResidueSeqenceCard* residues_sequence_;
            PdbResidueModificationCard* residue_modification_;
            PdbHeterogenCard* heterogens_;
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
