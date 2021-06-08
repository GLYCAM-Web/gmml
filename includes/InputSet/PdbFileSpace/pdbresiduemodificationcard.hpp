// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBRESIDUEMODIFICATIONCARD_HPP
#define PDBRESIDUEMODIFICATIONCARD_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbResidueModificationCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbResidueModificationCard();
            /*! \fn
              * Constructor with required parameters
              * @param id_code Identification code of a residue modification
              * @param residue_name Name of the residue in the residue modification
              * @param chain_id Chain identifier of the residue modification
              * @param sequence_number Sequence number of the residue modification
              * @param insertion_code Insertion code of the residue modification
              * @param standard_residue_name Standard name of the residue in the residue modification
              * @param dscr Short description for residue modification object
              */
            PdbResidueModificationCard(const std::string& id_code, const std::string& residue_name, char chain_id, int sequence_number,
                                   char insertion_code, const std::string& standard_residue_name, const std::string& dscr);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block Block of stream consists of residue modifications
              */
            PdbResidueModificationCard(std::stringstream& stream_block);
            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the id code in a residue modification
              * @return id_code_ attribute of the current object of this class
              */
            std::string GetIdCode();
            /*! \fn
              * An accessor function in order to access to the residue name in a residue modification
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the chain identifier in a residue modification
              * @return chain_id_ attribute of the current object of this class
              */
            char GetChainId();
            /*! \fn
              * An accessor function in order to access to the sequence number in a residue modification
              * @return sequence_number_ attribute of the current object of this class
              */
            int GetSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the insertion code in a residue modification
              * @return insertion_code_ attribute of the current object of this class
              */
            char GetInsertionCode();
            /*! \fn
              * An accessor function in order to access to the standard residue name in a residue modification
              * @return standard_residue_name_ attribute of the current object of this class
              */
            std::string GetStandardResidueName();
            /*! \fn
              * An accessor function in order to access to the description in a residue modification
              * @return dscr_ attribute of the current object of this class
              */
            std::string GetDscr();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the id code of the current object
              * Set the id_code_ attribute of the current residue modification
              * @param id_code The id code of the current object
              */
            void SetIdCode(const std::string id_code);
            /*! \fn
              * A mutator function in order to set the residue_name of the current object
              * Set the residue_name_ attribute of the current residue modification
              * @param residue_name The residue name of the current object
              */
            void SetResidueName(const std::string residue_name);
            /*! \fn
              * A mutator function in order to set the chain identifier of the current object
              * Set the chain_id_ attribute of the current residue modification
              * @param chain_id The chain identifier of the current object
              */
            void SetChainId(char chain_id);
            /*! \fn
              * A mutator function in order to set the sequence number of the current object
              * Set the sequence_number_ attribute of the current residue modification
              * @param sequence_number The sequence number of the current object
              */
            void SetSequenceNumber(int sequence_number);
            /*! \fn
              * A mutator function in order to set the insertion code of the current object
              * Set the insertion_code_ attribute of the current residue modification
              * @param insertion_code The insertion code of the current object
              */
            void SetInsertionCode(char insertion_code);
            /*! \fn
              * A mutator function in order to set the standard residue name of the current object
              * Set the standard_residue_name_ attribute of the current residue modification
              * @param standard_residue_name The standard residue name of the current object
              */
            void SetStandardResidueName(const std::string standard_residue_name);
            /*! \fn
              * A mutator function in order to set the describtion of the current object
              * Set the dscr_ attribute of the current residue modification
              * @param dscr The dscribtion of the current object
              */
            void SetDscr(const std::string dscr);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb residue modification contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string id_code_;                   /*!< Identification code of a residue modification >*/
            std::string residue_name_;              /*!< Name of the residue in a residue modification >*/
            char chain_id_;                         /*!< Chain identification of a residue in a residue modification >*/
            int sequence_number_;                   /*!< Sequence number of a residue in a residue modification >*/
            char insertion_code_;                   /*!< Insertion code of a residue in a residue modification >*/
            std::string standard_residue_name_;     /*!< Standard residue name of a residue in a residue modification >*/
            std::string dscr_;                      /*!< Short description for a residue in a residue modification >*/
    };
}

#endif // PDBRESIDUEMODIFICATIONCARD_HPP
