#ifndef PDBRESIDUESEQUENCE_HPP
#define PDBRESIDUESEQUENCE_HPP

#include <string>
#include <vector>

namespace PdbFileSpace
{
    class PdbResidueSequence
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbResidueSequence();
            PdbResidueSequence(char chain_id, int number_of_residues, const std::vector<std::string>& residue_names);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            char GetChainId();
            int GetNumberOfResidues();
            std::vector<std::string> GetResidueNames();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetChainId(char chain_id);
            void SetNumberOfResidues(int number_of_residues);
            void SetResidueNames(const std::vector<std::string> residue_names);
            void AddResidueName(const std::string residue_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            char chain_id_;
            int number_of_residues_;
            std::vector<std::string> residue_names_;
    };
}

#endif // PDBRESIDUESEQUENCE_HPP
