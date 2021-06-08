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
              * @param residue_name Name of site residue card record that appears in the first column of each line of a pdb file
              * @param residue_chain_id Chain identifier of a residue in a pdb file
              * @param residue_sequence_number Sequence number of a residue in a pdb file
              * @param residue_insertion_code Insertion code of a residue in a pdb file
              */
            PdbSiteResidue(const std::string &residue_name, char residue_chain_id, int residue_sequence_number, char residue_insertion_code);
            /*! \fn
              * Constructor with required parameters
              * @param section
              */
            PdbSiteResidue(const std::string& section);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
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
            int GetResidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the residue insertion code in a pdb site residue
              * @return resdiue_insertion_code_ attribute of the current object of this class
              */
            char GetResidueInsertionCode();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
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
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb site residue contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string residue_name_;          /*!< Name of the residue appears in site residue card of a pdb file >*/
            char residue_chain_id_;             /*!< Chain identifier of a residue that apears in site residue card of a pdb file >*/
            int residue_sequence_number_;       /*!< Sequence number of a residue that appears in a site residue card of a pdb file >*/
            char residue_insertion_code_;       /*!< Insertion code of a residue that appears in a site card residue card of a pdb file >*/
    };
}

#endif // PDBSITERESIDUE_HPP
