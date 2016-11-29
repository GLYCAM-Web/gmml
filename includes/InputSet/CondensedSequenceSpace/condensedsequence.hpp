#ifndef CONDENSEDSEQUENCE_HPP
#define CONDENSEDSEQUENCE_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <bitset>
#include "../../common.hpp"

namespace CondensedSequenceSpace
{
    class CondensedSequenceResidue;
    class CondensedSequenceAmberPrepResidue;
    // Options for condensed sequence: rotamers and glycosidic angles
    struct RotamersAndGlycosidicAnglesInfo{
        public:
            RotamersAndGlycosidicAnglesInfo(int linkage_index, std::vector<std::pair<std::string, std::vector<std::string> > > possible_rotamers,
                                            std::vector<std::pair<std::string, std::vector<std::string> > > selected_rotamers,
                                            std::vector<std::pair<std::string, double> > enabled_glycosidic_angles)
            {
                linkage_index_ = linkage_index;
                possible_rotamers_ = possible_rotamers;
                selected_rotamers_ = selected_rotamers;
                enabled_glycosidic_angles_ = enabled_glycosidic_angles;
            }

            int GetLinkageIndex(){
                return linkage_index_;
            }

            std::vector<std::pair<std::string, std::vector<std::string> > > GetPossibleRotamers(){
                return possible_rotamers_;
            }

            std::vector<std::pair<std::string, std::vector<std::string> > > GetSelectedRotamers(){
                return selected_rotamers_;
            }

            std::vector<std::pair<std::string, double> > GetEnabledGlycosidicAngles(){
                return enabled_glycosidic_angles_;
            }

            int linkage_index_;
            std::vector<std::pair<std::string, std::vector<std::string> > > possible_rotamers_;
            std::vector<std::pair<std::string, std::vector<std::string> > > selected_rotamers_;
            std::vector<std::pair<std::string, double> > enabled_glycosidic_angles_;
    };

    class CondensedSequence
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////         
            /*! \typedef
              * List of condensed sequence residues
              */
            typedef std::vector<CondensedSequenceResidue*> CondensedSequenceResidueVector;
            /*! \typedef
              * List of condensed sequence residues
              */
            typedef std::vector<gmml::CondensedSequenceTokenType> CondensedSequenceTokenTypeVector;
//            typedef std::vector<std::pair<CondensedSequenceResidue*, int> > CondensedSequenceResidueTree;
//            typedef std::vector<std::pair<CondensedSequenceAmberPrepResidue*, int> > CondensedSequenceAmberPrepResidueTree;

            typedef std::vector<CondensedSequenceResidue*> CondensedSequenceResidueTree;
            typedef std::vector<CondensedSequenceAmberPrepResidue*> CondensedSequenceAmberPrepResidueTree;       
            typedef std::pair<std::string, RotamersAndGlycosidicAnglesInfo*> RotamerNameInfoPair;
            typedef std::vector<RotamerNameInfoPair> CondensedSequenceRotamersAndGlycosidicAnglesInfo;
            typedef std::map<int, std::vector<std::vector<double> > > IndexLinkageConfigurationMap;
            typedef std::map<int, std::vector<std::vector<std::string> > > IndexConfigurationNameMap;
            typedef std::map<int, std::string> IndexNameMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            CondensedSequence();
            CondensedSequence(std::string sequence);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            CondensedSequenceResidueVector GetResidues();
            CondensedSequenceTokenTypeVector GetTokens();
            CondensedSequenceResidueTree GetCondensedSequenceResidueTree();
            CondensedSequenceAmberPrepResidueTree GetCondensedSequenceAmberPrepResidueTree();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetResidues(CondensedSequenceResidueVector residues);
            void AddResidue(CondensedSequenceResidue* residue);
            void SetTokens(CondensedSequenceTokenTypeVector tokens);
            void AddToken(gmml::CondensedSequenceTokenType token);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            int InsertNodeInCondensedSequenceResidueTree(CondensedSequenceResidue* condensed_residue, int parent_node_id = -1);
            int InsertNodeInCondensedSequenceAmberPrepResidueTree(CondensedSequenceAmberPrepResidue* condensed_amber_prep_residue, int parent_node_id = -1);
            void ParseCondensedSequence(std::string sequence);
            void BuildArrayTreeOfCondensedSequenceResidue();
            void BuildArrayTreeOfCondensedSequenceAmberPrepResidue(CondensedSequenceResidueTree residue_tree);
            std::string GetAmberPrepTerminalResidueCodeOfTerminalResidue(std::string terminal_residue_name);
            std::string GetAmberPrepResidueCodeOfCondensedResidue(CondensedSequenceResidue* condensed_residue, std::vector<int> open_valences, std::string parent_name);
            std::string GetFirstLetterOfAmberPrepResidueCode(std::bitset<10> open_valences_check);
            std::string GetSecondLetterOfAmberPrepResidueCode(std::string residue_name, std::string isomer);
            std::string GetThirdLetterOfAmberPrepResidueCode(std::string configuration, std::string ring_type);
            CondensedSequenceAmberPrepResidue* GetCondensedSequenceDerivativeAmberPrepResidue(std::string derivative_name, int derivative_index);
            CondensedSequenceRotamersAndGlycosidicAnglesInfo GetCondensedSequenceRotamersAndGlycosidicAnglesInfo(CondensedSequenceResidueTree residue_tree);
            int CountAllPossibleSelectedRotamers(CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info);
            int CountAllPossible28LinkagesRotamers(CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info);
            std::vector<std::vector<int> > CreateBaseMapAllPossibleSelectedRotamers(CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info);
            IndexLinkageConfigurationMap CreateIndexLinkageConfigurationMap(CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info,
                                                                            IndexNameMap& names);


            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the condensed sequence contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            CondensedSequenceResidueVector residues_;
            CondensedSequenceTokenTypeVector tokens_;
            CondensedSequenceResidueTree condensed_sequence_residue_tree_;
            CondensedSequenceAmberPrepResidueTree condensed_sequence_amber_prep_residue_tree_;

    };
}

#endif // CONDENSEDSEQUENCE_HPP
