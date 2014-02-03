// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

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
            /*! \fn
              * Default constructor
              */
            PdbSiteResidue();
            /*! \fn
              * Constructor with required parameters
              * @param residue_name
              * @param residue_chain_id
              * @param residue_sequence_number
              * @param residue_insertion_code
              */
            PdbSiteResidue(const std::string &residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code);
            PdbSiteResidue(const std::string& section);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue name in a pdb site residue
              * @return resdiue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the residue chain id in a pdb site residue
              * @return resdiue_chian_id_ attribute of the current object of this class
              */
            char GetResidueChainId();
            /*! \fn
              * An accessor function in order to access to the residue sequence number in a pdb site residue
              * @return resdiue_sequence_number_ attribute of the current object of this class
              */
            int GetresidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the residue insertion code in a pdb site residue
              * @return resdiue_insertion_code_ attribute of the current object of this class
              */
            char GetResidueInsertionCode();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current site residue
              * @param residue_name The residue name of the current object
              */
            void SetResidueName(const std::string residue_name);
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current site residue
              * @param residue_chain_id The residue chain id of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the residue sequence number of the current object
              * Set the residue_sequence_number_ attribute of the current site residue
              * @param residue_sequence_number The residue sequence number of the current object
              */
            void SetResidueSequenceNumber(int residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the residue insertion code of the current object
              * Set the residue_insertion_code_ attribute of the current site residue
              * @param residue_insertion_code The residue insertion code of the current object
              */
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
