#ifndef PDBSITERESIDUE_HPP
#define PDBSITERESIDUE_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbSiteResidue
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PdbSiteResidue();
            PdbSiteResidue(const std::string &residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code);
            PdbSiteResidue(const std::string& section);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetResidueName();
            char GetResidueChainId();
            int GetresidueSequenceNumber();
            char GetResidueInsertionCode();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetResidueName(const std::string residue_name);
            void SetResidueChainId(char residue_chain_id);
            void SetResidueSequenceNumber(int residue_sequence_number);
            void SetResidueInsertionCode(char residue_insertion_code);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string residue_name_;
            char residue_chain_id_;
            int residue_sequence_number_;
            char residue_insertion_code_;
    };
}

#endif // PDBSITERESIDUE_HPP
